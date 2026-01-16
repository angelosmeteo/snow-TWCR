
# ============================================================
# snow-TWCR | run_all.R  (MAIN ENTRY POINT)
#
# This script SOURCES helper files from /code and runs the full workflow.
# It automatically uses all AWS stations that have valid lat/lon/elev_m.
#
# REQUIREMENTS (user provides locally):
#   - Excel file: snow_field_measurements.xlsx
#     (sheets: depth_raw, density, stations_meta, stations_10min)
#   - DEM GeoTIFF: modi.tif (must have CRS)
#
# OUTPUTS:
#   - maps_final/depth/Depth_YYYYMMDD_raw.tif
#   - maps_final/depth/Depth_YYYYMMDD_capXX.tif
#   - maps_final/swe/SWE_YYYYMMDD_capXX.tif
#   - maps_final/met/hillshade_UTM.tif (optional)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(sf)
  library(terra)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(lubridate)
  library(tidyr)
})

# -----------------------------
# SOURCE MODULES (MUST EXIST)
# -----------------------------
source(file.path("code", "helpers_io.R"))
source(file.path("code", "aws_psnow.R"))
source(file.path("code", "mapping_model.R"))

# -----------------------------
# USER PATHS (EDIT THESE)
# -----------------------------
base_dir  <- "C:/Users/AngelosT/Documents/DSS"
xlsx_path <- file.path(base_dir, "snow_field_measurements.xlsx")
dem_path  <- file.path(base_dir, "modi.tif")

# One clean output folder
out_base  <- file.path(base_dir, "maps_final")
out_depth <- file.path(out_base, "depth")
out_swe   <- file.path(out_base, "swe")
out_met   <- file.path(out_base, "met")

dir.create(out_depth, showWarnings = FALSE, recursive = TRUE)
dir.create(out_swe,   showWarnings = FALSE, recursive = TRUE)
dir.create(out_met,   showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# STUDY REGION (EPSG:4326 bbox)
# -----------------------------
bbox_ll <- list(S=40.612980, N=40.693127, W=23.600922, E=23.747495)

# UTM zone 34N
utm_crs <- "EPSG:32634"

# -----------------------------
# FIXED DEPTH CAPS (REQUESTED)
# -----------------------------
DEPTH_CAP_FIXED <- c("2025-01-14" = 80,
                     "2025-01-15" = 150)

# -----------------------------
# MODEL SETTINGS
# -----------------------------
P0_CUTOFF <- 0.06
P_FLOOR_HIGH <- 0.85
E_MARGIN <- 10

# Storm + AWS settings
TZ_USE <- "Europe/Athens"
STORM_START_LOCAL <- "2025-01-12 00:00:00"
SURVEY_TIME_LOCAL <- "12:00:00"

LAPSE_C_PER_KM <- -6.5
IDW_POWER <- 2

# Wet-bulb transition band centered near rain–snow transition
TW_CENTER <- 0.5
TW_SNOW_FULL <- TW_CENTER - 1.0   # -0.5°C
TW_RAIN_FULL <- TW_CENTER + 1.0   #  1.5°C

# Dewpoint QC
DEWPOINT_QC_THRESH <- 2.0
APPLY_DEWPOINT_QC  <- TRUE

# Optional hillshade output (handy in QGIS)
WRITE_HILLSHADE <- TRUE

# ============================================================
# 1) READ FIELD MEASUREMENTS
# ============================================================
depth_raw <- read_excel(xlsx_path, sheet="depth_raw") %>%
  mutate(survey_date = as.Date(survey_date))

density_raw <- read_excel(xlsx_path, sheet="density") %>%
  mutate(survey_date = as.Date(survey_date))

# Point-level depth mean from replicates
depth_pt <- depth_raw %>%
  group_by(survey_date, point_id, lat, lon, zone) %>%
  summarise(depth_mean_cm = mean(depth_cm, na.rm=TRUE),
            depth_sd_cm   = sd(depth_cm, na.rm=TRUE),
            n_depth       = dplyr::n(),
            .groups="drop")

# ============================================================
# 2) DEM PREP: crop -> UTM
# ============================================================
dem <- rast(dem_path)
if (is.na(crs(dem))) stop("❌ DEM has no CRS. Assign CRS and re-save.")

bbox_poly_ll <- st_as_sfc(st_bbox(c(xmin=bbox_ll$W, ymin=bbox_ll$S,
                                    xmax=bbox_ll$E, ymax=bbox_ll$N), crs=4326))

dem_ll <- mask(crop(dem, vect(bbox_poly_ll), snap="out"), vect(bbox_poly_ll))
if (ncell(dem_ll) == 0) stop("❌ DEM crop produced 0 cells. Check bbox overlap.")

dem_c <- project(dem_ll, utm_crs, method="bilinear")
names(dem_c) <- "elev"

if (WRITE_HILLSHADE) {
  hs <- shade(terrain(dem_c, "slope", unit="radians"),
              terrain(dem_c, "aspect", unit="radians"),
              angle=45, direction=315)
  writeRaster(hs, file.path(out_met, "hillshade_UTM.tif"), overwrite=TRUE)
}

# ============================================================
# 3) DISTANCE TO SEA (Natural Earth land polygons)
# ============================================================
dist_sea_km <- build_distance_to_sea_km(dem_c, bbox_poly_ll, utm_crs)
names(dist_sea_km) <- "dist_sea_km"

# ============================================================
# 4) GRID DF (DEM cells only)
# ============================================================
grid_df <- as.data.frame(dem_c, xy=TRUE, cells=TRUE, na.rm=TRUE)
idx <- grid_df$cell

grid_df$dist_sea_km <- values(dist_sea_km)[idx]
grid_df$dist_sea_km[is.na(grid_df$dist_sea_km)] <- median(grid_df$dist_sea_km, na.rm=TRUE)
elev_vec <- grid_df$elev

# ============================================================
# 5) BUILD SWE POINT DATA (density imputation by elevation)
# ============================================================
# robust to missing coords: returns NA elev for those points instead of crashing
depth_pt$elev_dem <- extract_elev_points(depth_pt %>% select(lon, lat), dem_c)
density_raw$elev_dem <- extract_elev_points(density_raw %>% select(lon, lat), dem_c)

dens_global_mean <- mean(density_raw$density_gcm3, na.rm=TRUE)

dat_swe <- build_swe_points_from_depth_density(
  depth_pt = depth_pt,
  density_raw = density_raw,
  dem_elev_field = "elev_dem",
  dens_global_mean = dens_global_mean
)

# ============================================================
# 6) READ + QC AWS AND PRECOMPUTE IDW WEIGHTS (AUTOMATIC ALL STATIONS)
# ============================================================
aws <- read_and_prepare_aws(
  xlsx_path = xlsx_path,
  TZ_USE = TZ_USE,
  DEWPOINT_QC_THRESH = DEWPOINT_QC_THRESH,
  APPLY_DEWPOINT_QC = APPLY_DEWPOINT_QC
)

stations_meta  <- aws$stations_meta
stations_10min <- aws$stations_10min

cat("\n✅ Stations used (valid lat/lon/elev_m):\n")
print(stations_meta$station_id)

Wnorm <- build_station_weights_Wnorm(
  dem_c = dem_c,
  grid_df = grid_df,
  stations_meta = stations_meta,
  utm_crs = utm_crs,
  IDW_POWER = IDW_POWER
)

# ============================================================
# 7) RUN DAYS (Depth + SWE)
# ============================================================
run_one_day <- function(day_date) {

  dtag <- format(as.Date(day_date), "%Y-%m-%d")
  if (!(dtag %in% names(DEPTH_CAP_FIXED))) stop(paste0("❌ Missing fixed cap for ", dtag))
  depth_cap_day <- as.numeric(DEPTH_CAP_FIXED[[dtag]])

  # storm-integrated Psnow raster for this day
  Psnow_r <- build_storm_psnow_surface(
    day_date = as.Date(day_date),
    stations_10min = stations_10min,
    stations_meta  = stations_meta,
    Wnorm = Wnorm,
    dem_c = dem_c,
    grid_df = grid_df,
    TZ_USE = TZ_USE,
    STORM_START_LOCAL = STORM_START_LOCAL,
    SURVEY_TIME_LOCAL = SURVEY_TIME_LOCAL,
    LAPSE_C_PER_KM = LAPSE_C_PER_KM,
    TW_SNOW_FULL = TW_SNOW_FULL,
    TW_RAIN_FULL = TW_RAIN_FULL
  )

  # ---- DEPTH ----
  depth_pts_day <- depth_pt %>%
    filter(survey_date == as.Date(day_date)) %>%
    mutate(val = depth_mean_cm)

  # raw depth
  predict_two_stage_to_tif(
    points_df = depth_pts_day,
    varname = "Depth",
    dem_c = dem_c,
    dist_sea_km = dist_sea_km,
    Psnow_r = Psnow_r,
    grid_df = grid_df,
    P0_CUTOFF = P0_CUTOFF,
    P_FLOOR_HIGH = P_FLOOR_HIGH,
    E_MARGIN = E_MARGIN,
    cap_value = NA,
    out_tif = file.path(out_depth, paste0("Depth_", format(as.Date(day_date), "%Y%m%d"), "_raw.tif"))
  )

  # capped depth
  predict_two_stage_to_tif(
    points_df = depth_pts_day,
    varname = "Depth",
    dem_c = dem_c,
    dist_sea_km = dist_sea_km,
    Psnow_r = Psnow_r,
    grid_df = grid_df,
    P0_CUTOFF = P0_CUTOFF,
    P_FLOOR_HIGH = P_FLOOR_HIGH,
    E_MARGIN = E_MARGIN,
    cap_value = depth_cap_day,
    out_tif = file.path(out_depth, paste0("Depth_", format(as.Date(day_date), "%Y%m%d"), "_cap", round(depth_cap_day), ".tif"))
  )

  # ---- SWE ----
  swe_pts_day <- dat_swe %>%
    filter(survey_date == as.Date(day_date)) %>%
    mutate(val = swe_mm)

  dens_day_vec <- density_raw %>% filter(survey_date == as.Date(day_date)) %>% pull(density_gcm3)
  dens_ref <- if (length(dens_day_vec) >= 3) as.numeric(quantile(dens_day_vec, 0.95, na.rm=TRUE)) else dens_global_mean
  SWE_CAP_DAY <- depth_cap_day * dens_ref * 10

  predict_two_stage_to_tif(
    points_df = swe_pts_day,
    varname = "SWE",
    dem_c = dem_c,
    dist_sea_km = dist_sea_km,
    Psnow_r = Psnow_r,
    grid_df = grid_df,
    P0_CUTOFF = P0_CUTOFF,
    P_FLOOR_HIGH = P_FLOOR_HIGH,
    E_MARGIN = E_MARGIN,
    cap_value = SWE_CAP_DAY,
    out_tif = file.path(out_swe, paste0("SWE_", format(as.Date(day_date), "%Y%m%d"), "_cap", round(SWE_CAP_DAY), ".tif"))
  )

  cat("✅ Finished day:", as.character(day_date), "\n")
}

run_one_day(as.Date("2025-01-14"))
run_one_day(as.Date("2025-01-15"))

cat("\n✅ Done.\nAll outputs under:", out_base, "\n")

