
# snow-TWCR

**snow-TWCR** is an R workflow for mapping **snow depth (cm)** and **snow water equivalent (SWE, mm)** across a study area using (i) field snow observations and (ii) storm-integrated meteorological forcing from automatic weather stations (AWS). The workflow produces GeoTIFF rasters on a DEM grid (ready for QGIS).

**Repository title (for citation):**  
*snow-TWCR: Storm-integrated mapping of snow depth and SWE using wet-bulb temperature and a two-stage occurrence–magnitude model (field-conditioned)*

---

## What this repository includes

- **`code/`**  
  R scripts that:
  - compute storm-cumulative **snow‑favorable precipitation** (`P_snow`) from AWS precipitation and wet‑bulb temperature (Tw), adjusted by a constant atmospheric lapse rate and blended across stations using distance weighting
  - map **snow depth** and **SWE** over the full study area using a two-stage occurrence × magnitude model
  - export final GeoTIFFs to a single clean output folder (`maps_final/`)

- **`formulas/`**  
  - `model_formulas.pdf` containing the full model equations and definitions used in the workflow.

- **`CITATION.cff`**  
  Citation metadata so the repository can be cited correctly (GitHub shows “Cite this repository”).

---

## What this repository does NOT include (by design)

This repository intentionally does **not** include:
- the raw field observation files
- AWS raw data files
- the DEM raster(s)

The workflow uses field measurements as **inputs**, but the measurement datasets are not distributed in this repository. This repo provides a **code + formulas** package accompanying a manuscript.

---

## Method summary (plain language)

### 1) Storm-integrated snow-favorable precipitation (`P_snow`)
From the storm start (**12 Jan 2025 00:00 local**) to each survey time (**12:00 local**) at **10-minute resolution**, I estimate at every DEM grid cell:
- precipitation as a **distance-weighted blend** of all available AWS stations (continuous weighting; no nearest-station boundary)
- wet-bulb temperature (Tw) from air temperature and relative humidity, adjusted to elevation using a **constant lapse rate** (−6.5 °C/km)
- a **snow fraction** from Tw using a transition band (snow → mixed → rain)

I then accumulate snow-favorable precipitation:
\[
P_{\text{snow}}(x) = \sum_t P(x,t)\,f_s(x,t)
\]

### 2) Two-stage mapping (occurrence × magnitude)
I map snow depth and SWE using a two-stage model:
1. **Occurrence**: logistic regression predicting snow presence/absence  
2. **Magnitude**: linear regression on `log(1 + value)` fitted only where snow exists, then back-transformed  
Final prediction is **probability × magnitude**, with additional rules to remove spurious low-elevation snow, enforce realistic high-elevation snow probability, and hard-condition predictions at observation grid cells (field-conditioned).

---

## Weather stations (automatic selection)

The workflow automatically uses **all AWS stations** that have valid metadata in `stations_meta` (**latitude, longitude, and elevation**). Stations missing any of these fields are ignored (because distances and lapse-rate adjustments cannot be computed).

---

## Requirements

- R (≥ 4.x recommended)

R packages used:
- `readxl`, `dplyr`, `sf`, `terra`, `lubridate`, `tidyr`, `rnaturalearth`, `rnaturalearthdata`
- Optional: `brglm2` (used if installed for bias-reduced logistic regression when logistic separation occurs)

---

## How to run (if you provide input data locally)

1. Open `code/run_all.R`
2. Edit the local paths:
   - `base_dir`
   - `xlsx_path` (Excel with field + AWS sheets)
   - `dem_path`  (DEM GeoTIFF with valid CRS)
3. Run:
   - in R: `source("code/run_all.R")`
   - or command line: `Rscript code/run_all.R`

---

## Outputs

The workflow writes GeoTIFF rasters to:

- `maps_final/depth/`
  - `Depth_YYYYMMDD_raw.tif`
  - `Depth_YYYYMMDD_capXX.tif`

- `maps_final/swe/`
  - `SWE_YYYYMMDD_capXX.tif`

- `maps_final/met/`
  - `hillshade_UTM.tif` (optional background for QGIS)

For the January 2025 case study, depth caps were fixed to:
- 2025‑01‑14: **80 cm**
- 2025‑01‑15: **150 cm**

SWE caps are derived from the fixed depth caps and a robust density estimate (see `model_formulas.pdf`).

---

## Citation / DOI

- If you use or adapt this code, please cite it using the information in **`CITATION.cff`**.

Placeholders:
- Zenodo Version DOI (exact release used in the paper): 10.5281/zenodo.18270203
- Zenodo Concept DOI (all versions): 10.5281/zenodo.18270203

Repository URL: https://github.com/angelosmeteo/snow-TWCR

---

## License

See the `LICENSE` file.

---

## Contact

Angelos Theodorou

https://theodorouweather.com/ 
