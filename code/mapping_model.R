# ============================================================
# mapping_model.R
# Two-stage mapping:
#   1) logistic regression for snow presence
#   2) linear regression on log1p(value) for magnitude (snow-present only)
# Soft coupling: yhat = p_hat * m_hat
# Post-processing:
#   - low-elev cutoff below snowline proxy
#   - high-elev probability floor
# Hard conditioning at measurement cells
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
})

predict_two_stage_to_tif <- function(points_df, varname,
                                     dem_c, dist_sea_km, Psnow_r, grid_df,
                                     P0_CUTOFF, P_FLOOR_HIGH, E_MARGIN,
                                     cap_value = NA,
                                     out_tif) {

  idx <- grid_df$cell

  # attach Psnow predictor to grid_df
  grid_df$Psnow_storm_mm <- values(Psnow_r)[idx]
  grid_df$Psnow_storm_mm[is.na(grid_df$Psnow_storm_mm)] <- median(grid_df$Psnow_storm_mm, na.rm=TRUE)

  # UTM points
  utm_crs <- crs(dem_c)
  pts_sf <- st_transform(st_as_sf(points_df, coords=c("lon","lat"), crs=4326), utm_crs)
  pts_v  <- vect(pts_sf)
  xy <- crds(pts_v)

  elev_pts <- terra::extract(dem_c, pts_v)[,2]
  dsea_pts <- terra::extract(dist_sea_km, pts_v)[,2]
  psno_pts <- terra::extract(Psnow_r, pts_v)[,2]

  pts <- data.frame(points_df,
                    x=xy[,1], y=xy[,2],
                    elev=elev_pts,
                    dist_sea_km=dsea_pts,
                    Psnow_storm_mm=psno_pts) %>%
    filter(!is.na(elev))

  pts$snow_present <- ifelse(pts$val > 0, 1, 0)

  # Stage 1: presence model
  pres_formula <- snow_present ~ elev + I(elev^2) + log1p(dist_sea_km) + log1p(Psnow_storm_mm)

  if (requireNamespace("brglm2", quietly=TRUE)) {
    pres_fit <- brglm2::brglm(pres_formula, data=pts, family=binomial())
  } else {
    pres_fit <- glm(pres_formula, data=pts, family=binomial(), control=list(maxit=100))
  }

  p_snow <- predict(pres_fit, newdata=grid_df, type="response")

  # snowline proxy
  e_snow <- pts$elev[pts$snow_present == 1]
  E_snowline <- if (length(e_snow) >= 3) as.numeric(quantile(e_snow, 0.10, na.rm=TRUE)) else 200

  # high-elevation probability floor
  if (length(e_snow) >= 3) {
    E_floor <- as.numeric(quantile(e_snow, 0.10, na.rm=TRUE))
    high_idx <- which(grid_df$elev >= (E_floor + E_MARGIN))
    p_snow[high_idx] <- pmax(p_snow[high_idx], P_FLOOR_HIGH)
  }

  # Stage 2: magnitude model (snow-present only)
  pts_pos <- pts %>% filter(snow_present == 1)
  mag_formula <- log1p(val) ~ elev + I(elev^2) + log1p(dist_sea_km) + log1p(Psnow_storm_mm)
  mag_fit <- lm(mag_formula, data=pts_pos, na.action=na.exclude)

  mag_pred <- pmax(0, expm1(predict(mag_fit, newdata=grid_df)))

  # Soft coupling
  pred <- p_snow * mag_pred

  # Low-elevation cutoff below snowline proxy
  low_idx  <- which(grid_df$elev < E_snowline)
  kill_idx <- low_idx[p_snow[low_idx] < P0_CUTOFF]
  pred[kill_idx] <- 0
  pred <- pmax(0, pred)

  # Apply cap if provided
  if (!is.na(cap_value)) pred <- pmin(pred, cap_value)

  # Rasterize
  r <- rast(dem_c); values(r) <- NA_real_
  r[idx] <- pred

  # Hard conditioning at measurement cells
  pt_cells <- cellFromXY(dem_c, cbind(pts$x, pts$y))
  obs_cell <- data.frame(cell=pt_cells, obs=pts$val)
  obs_cell <- obs_cell[!is.na(obs_cell$cell), ]
  obs_cell <- aggregate(obs ~ cell, data=obs_cell, FUN=mean)
  r[obs_cell$cell] <- obs_cell$obs

  names(r) <- varname
  writeRaster(r, out_tif, overwrite=TRUE)
  cat("✅ Wrote:", out_tif, "\n")

  invisible(r)
}

