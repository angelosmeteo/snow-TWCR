# ============================================================
# aws_psnow.R
# AWS processing:
#  - Automatically uses ALL stations with valid lat/lon/elev_m from stations_meta
#  - Filters stations_10min to only those valid stations
#  - Computes wet-bulb Tw from T and RH
#  - Optional dewpoint QC (Td from T/RH vs observed Td)
# Storm Psnow:
#  - distance-weighted precip and Tw0 across ALL valid stations
#  - lapse-rate adjustment to each DEM cell elevation
#  - snow fraction from Tw transition band
#  - accumulate Psnow over storm period (10-min)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(sf)
  library(terra)
  library(lubridate)
  library(tidyr)
})

# ---- dewpoint from T & RH (Magnus)
dewpoint_from_T_RH <- function(T_C, RH_pct) {
  a <- 17.27; b <- 237.7
  RH <- pmax(pmin(RH_pct, 100), 1) / 100
  gamma <- log(RH) + (a * T_C) / (b + T_C)
  (b * gamma) / (a - gamma)
}

# ---- wet-bulb from T & RH (Stull-type approximation)
wetbulb_from_T_RH <- function(T_C, RH_pct) {
  RH <- pmax(pmin(RH_pct, 100), 1)
  T_C * atan(0.151977 * sqrt(RH + 8.313659)) +
    atan(T_C + RH) - atan(RH - 1.676331) +
    0.00391838 * RH^(3/2) * atan(0.023101 * RH) - 4.686035
}

# ---- Excel time parser (handles 1900 times)
parse_time_to_hms <- function(x, TZ_USE) {
  if (is.na(x)) return("00:00:00")
  s <- as.character(x)
  if (s %in% c("0", "0.0")) return("00:00:00")
  if (grepl("1900", s)) {
    tt <- suppressWarnings(as.POSIXct(s, tz=TZ_USE))
    if (!is.na(tt)) return(format(tt, "%H:%M:%S"))
  }
  if (grepl("^\\d{1,2}:\\d{2}(:\\d{2})?$", s)) {
    if (nchar(s) == 5) return(paste0(s, ":00"))
    return(s)
  }
  tt <- suppressWarnings(as.POSIXct(s, tz=TZ_USE))
  if (!is.na(tt)) return(format(tt, "%H:%M:%S"))
  "00:00:00"
}

# ---- read and QC AWS (AUTOMATIC all valid stations)
read_and_prepare_aws <- function(xlsx_path, TZ_USE, DEWPOINT_QC_THRESH=2, APPLY_DEWPOINT_QC=TRUE) {

  # 1) Read ALL stations that have valid metadata (lat/lon/elev_m)
  stations_meta <- readxl::read_excel(xlsx_path, sheet="stations_meta") %>%
    mutate(elev_m = as.numeric(elev_m)) %>%
    filter(!is.na(lat), !is.na(lon), !is.na(elev_m))

  if (nrow(stations_meta) < 1) {
    stop("❌ No valid stations in stations_meta (need lat, lon, elev_m).")
  }

  # 2) Read 10-min data and keep only those station IDs
  stations_10min <- readxl::read_excel(xlsx_path, sheet="stations_10min") %>%
    mutate(
      date = as.Date(date),
      time_str = vapply(time, parse_time_to_hms, character(1), TZ_USE=TZ_USE),
      datetime = as.POSIXct(paste(date, time_str), tz=TZ_USE),
      temp_C = as.numeric(temp_C),
      rh_pct = as.numeric(rh_pct),
      dewpoint_C = as.numeric(dewpoint_C),
      precip_mm = as.numeric(precip_mm)
    ) %>%
    filter(!is.na(datetime)) %>%
    filter(station_id %in% stations_meta$station_id) %>%
    mutate(
      dew_calc_C = dewpoint_from_T_RH(temp_C, rh_pct),
      dew_absdiff = abs(dew_calc_C - dewpoint_C),
      Tw_C = wetbulb_from_T_RH(temp_C, rh_pct)
    )

  if (APPLY_DEWPOINT_QC) {
    stations_10min <- stations_10min %>%
      filter(is.na(dew_absdiff) | dew_absdiff <= DEWPOINT_QC_THRESH)
  }

  list(stations_meta=stations_meta, stations_10min=stations_10min)
}

# ---- build normalized inverse-distance weights Wnorm (cells x stations)
build_station_weights_Wnorm <- function(dem_c, grid_df, stations_meta, utm_crs, IDW_POWER=2) {

  idx <- grid_df$cell

  st_sf <- st_transform(st_as_sf(stations_meta, coords=c("lon","lat"), crs=4326), utm_crs)
  st_v  <- vect(st_sf)

  nst <- nrow(stations_meta)
  W <- matrix(NA_real_, nrow=length(idx), ncol=nst)

  for (i in 1:nst) {
    d_km <- terra::distance(dem_c, st_v[i]) / 1000
    dval <- terra::values(d_km)[idx]
    dval[dval == 0] <- 0.001
    W[,i] <- 1 / (dval ^ IDW_POWER)
  }

  Wsum <- rowSums(W, na.rm=TRUE)
  sweep(W, 1, Wsum, FUN="/")
}

# ---- snow fraction from wet-bulb (linear transition)
snowfrac_from_Tw <- function(Tw, TW_SNOW_FULL, TW_RAIN_FULL) {
  out <- rep(NA_real_, length(Tw))
  out[Tw <= TW_SNOW_FULL] <- 1
  out[Tw >= TW_RAIN_FULL] <- 0
  mid <- which(Tw > TW_SNOW_FULL & Tw < TW_RAIN_FULL)
  out[mid] <- 1 - (Tw[mid] - TW_SNOW_FULL) / (TW_RAIN_FULL - TW_SNOW_FULL)
  pmax(0, pmin(1, out))
}

# ---- storm-integrated Psnow raster for a given day
build_storm_psnow_surface <- function(day_date,
                                      stations_10min, stations_meta, Wnorm,
                                      dem_c, grid_df,
                                      TZ_USE, STORM_START_LOCAL, SURVEY_TIME_LOCAL,
                                      LAPSE_C_PER_KM,
                                      TW_SNOW_FULL, TW_RAIN_FULL) {

  idx <- grid_df$cell
  elev_vec <- grid_df$elev

  storm_start <- as.POSIXct(STORM_START_LOCAL, tz=TZ_USE)
  survey_dt   <- as.POSIXct(paste(as.Date(day_date), SURVEY_TIME_LOCAL), tz=TZ_USE)

  sub <- stations_10min %>%
    filter(datetime >= storm_start, datetime <= survey_dt) %>%
    select(datetime, station_id, precip_mm, Tw_C)

  if (nrow(sub) == 0) stop("❌ No AWS data found in storm window for this day.")

  st_ids <- stations_meta$station_id

  # full time x station grid, fill missing precip with 0, fill Tw by down/up
  all_times <- sort(unique(sub$datetime))
  sub_full <- tidyr::expand_grid(datetime=all_times, station_id=st_ids) %>%
    left_join(sub, by=c("datetime","station_id")) %>%
    arrange(station_id, datetime) %>%
    group_by(station_id) %>%
    tidyr::fill(Tw_C, .direction="downup") %>%
    ungroup() %>%
    mutate(precip_mm = ifelse(is.na(precip_mm), 0, precip_mm))

  precip_w <- sub_full %>%
    select(datetime, station_id, precip_mm) %>%
    pivot_wider(names_from=station_id, values_from=precip_mm, values_fill=0) %>%
    arrange(datetime)

  Tw_w <- sub_full %>%
    select(datetime, station_id, Tw_C) %>%
    pivot_wider(names_from=station_id, values_from=Tw_C) %>%
    arrange(datetime)

  precip_mat <- as.matrix(precip_w[, st_ids, drop=FALSE])
  Tw_mat     <- as.matrix(Tw_w[, st_ids, drop=FALSE])

  z_st <- stations_meta$elev_m

  # convert station Tw to sea-level equivalent Tw0
  Tw0_mat <- sweep(Tw_mat, 2, (LAPSE_C_PER_KM * z_st/1000), FUN="-")

  Psnow <- rep(0, length(idx))

  for (k in 1:nrow(precip_mat)) {

    p_k   <- precip_mat[k, ]
    tw0_k <- Tw0_mat[k, ]

    # distance-weighted precip at each cell
    p_grid <- as.numeric(Wnorm %*% p_k)

    # distance-weighted Tw0 then lapse to cell elevation
    tw0_grid <- as.numeric(Wnorm %*% tw0_k)
    tw_grid  <- tw0_grid + LAPSE_C_PER_KM * (elev_vec/1000)

    sf_grid <- snowfrac_from_Tw(tw_grid, TW_SNOW_FULL, TW_RAIN_FULL)
    Psnow <- Psnow + (p_grid * sf_grid)
  }

  rPsnow <- rast(dem_c); values(rPsnow) <- NA_real_
  rPsnow[idx] <- Psnow
  names(rPsnow) <- "Psnow_storm_mm"
  rPsnow
}

