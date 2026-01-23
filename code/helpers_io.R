
# ============================================================
# helpers_io.R
# Utilities:
#   - extract DEM elevation at points (robust to missing coords)
#   - distance-to-sea raster (Natural Earth PHYSICAL coastline; excludes lakes)
#   - SWE point builder with elevation-based density interpolation
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

# ---- extract DEM elevation at lon/lat points (robust to NA coords)
extract_elev_points <- function(df_lonlat, dem_rast) {

  if (!all(c("lon","lat") %in% names(df_lonlat))) {
    stop("df_lonlat must contain columns: lon, lat")
  }

  out <- rep(NA_real_, nrow(df_lonlat))
  ok <- is.finite(df_lonlat$lon) & is.finite(df_lonlat$lat)
  if (!any(ok)) return(out)

  pts_sf <- sf::st_as_sf(df_lonlat[ok, ], coords=c("lon","lat"), crs=4326, remove=FALSE)
  pts_sf <- sf::st_transform(pts_sf, terra::crs(dem_rast))
  pts_v  <- terra::vect(pts_sf)

  out_ok <- terra::extract(dem_rast, pts_v)[,2]
  out[ok] <- out_ok
  out
}

# ---- distance-to-sea raster (km) using Natural Earth PHYSICAL coastline
#      This avoids inland lakes (e.g., Lake Volvi) being treated as "sea".
#
# Notes:
# - Uses Natural Earth "coastline" (physical) as the target geometry.
# - Distance is computed in projected units (UTM) and returned in km.
# - terra::distance() measures distance to the supplied vector features. [2](https://github.com/topics/snow)
# - Natural Earth provides coastline separately from lakes/reservoirs. [1](https://tteck.github.io/Proxmox/)
build_distance_to_sea_km <- function(dem_c, bbox_poly_ll, utm_crs) {

  sf::sf_use_s2(FALSE)

  # Get marine coastline (physical dataset; lakes are separate theme)
  coast_sf <- rnaturalearth::ne_download(
    scale = 50, category = "physical", type = "coastline", returnclass = "sf"
  )

  # Ensure CRS compatibility for intersection
  if (sf::st_crs(coast_sf) != sf::st_crs(bbox_poly_ll)) {
    coast_sf <- sf::st_transform(coast_sf, sf::st_crs(bbox_poly_ll))
  }

  # Clip coastline to bbox
  coast_clip <- suppressWarnings(sf::st_intersection(coast_sf, bbox_poly_ll))

  if (nrow(coast_clip) == 0) {
    stop("❌ coastline clip empty. Check bbox extent/CRS.")
  }

  # Make valid and union for robustness (optional but safe)
  coast_union <- tryCatch(
    sf::st_union(sf::st_make_valid(coast_clip)),
    error = function(e) sf::st_union(coast_clip)
  )

  # Project to UTM and compute Euclidean distance (m -> km)
  coast_v <- terra::project(terra::vect(coast_union), utm_crs)

  terra::distance(dem_c, coast_v) / 1000
}

# ---- monotonic decreasing density with elevation (isotonic regression)
estimate_density_by_elev <- function(elev_query, elev_obs, dens_obs) {

  ok <- is.finite(elev_obs) & is.finite(dens_obs)
  elev_obs <- elev_obs[ok]
  dens_obs <- dens_obs[ok]

  if (length(dens_obs) == 0) return(rep(NA_real_, length(elev_query)))
  if (length(dens_obs) == 1) return(rep(dens_obs[1], length(elev_query)))

  ord <- order(elev_obs)
  x <- elev_obs[ord]
  y <- dens_obs[ord]

  # isoreg enforces increasing; apply to (-y) to enforce decreasing in y
  iso <- stats::isoreg(x, -y)
  y_mono <- -iso$yf

  # interpolate and clamp beyond observed range (rule=2)
  y_hat <- approx(x=x, y=y_mono, xout=elev_query, rule=2)$y

  # clamp to fitted min/max
  y_hat <- pmax(min(y_mono, na.rm=TRUE), pmin(max(y_mono, na.rm=TRUE), y_hat))
  y_hat
}

# ---- build SWE points from depth points + (partial) density measurements
build_swe_points_from_depth_density <- function(depth_pt, density_raw, dem_elev_field="elev_dem", dens_global_mean) {

  out <- depth_pt %>%
    left_join(density_raw %>% dplyr::select(survey_date, point_id, density_gcm3),
              by=c("survey_date","point_id")) %>%
    group_by(survey_date) %>%
    group_modify(function(d, key) {

      day <- key$survey_date[1]

      dens_day <- density_raw %>%
        filter(survey_date == day) %>%
        filter(!is.na(density_gcm3), !is.na(.data[[dem_elev_field]]))

      if (nrow(dens_day) >= 1) {
        elev_obs <- dens_day[[dem_elev_field]]
        dens_obs <- dens_day$density_gcm3
        dens_est <- estimate_density_by_elev(d[[dem_elev_field]], elev_obs, dens_obs)
      } else {
        dens_est <- rep(dens_global_mean, nrow(d))
      }

      d$density_used <- d$density_gcm3
      d$density_used[is.na(d$density_used)] <- dens_est[is.na(d$density_used)]

      # if depth <= 0, density and SWE should be 0
      d$density_used <- ifelse(d$depth_mean_cm <= 0, 0, d$density_used)
      d$swe_mm <- ifelse(d$depth_mean_cm <= 0, 0, d$depth_mean_cm * d$density_used * 10)

      d
    }) %>%
    ungroup()

  out
}
