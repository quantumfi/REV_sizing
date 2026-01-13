# Sizing the Representative Element Volume of Cylindrical CT Scans  
## Supplementary Material: Software Documentation

**Authors:** Fernando Alonso-Marroquin, Abdullah Alqubalee, and Christian Tantardini

---

## Software purpose and scope

This MATLAB software is intended for the analysis of long cylindrical core samples (e.g., μCT volumes or other 3D imaging modalities) in which spatial properties such as porosity, phase fraction, or burrow-related indicators may exhibit slow axial variation.

In sufficiently long specimens, these low-frequency changes can act as a non-stationary trend along the core axis and can bias statistical descriptors that implicitly assume stationarity. The software provides a reproducible workflow to:

1. **Preprocess cylindrical 3D datasets** (binarization and enforcement of cylindrical support),
2. **Separate slow axial trends** from residual fluctuations using rolling-window detrending,
3. **Quantify whether the residual field** is compatible with an approximately stationary description over the scales of interest.

Two complementary analyses are implemented:

- **Axial analysis (1D reduction of 3D volume):**  
  The 3D volume is reduced to a slice-wise signal (e.g., mean phase fraction per slice) and a moving-average detrending is applied over a range of window sizes. Candidate windows are evaluated using the **excess kurtosis** of the detrended residual as a compact indicator of Gaussian-likeness, supporting a practical choice of detrending scale.

- **Transverse (slice-wise) spatial analysis (2D within a slice):**  
  Representativeness within a single cross-sectional slice is evaluated by computing **disk-averaged** estimates of the target property as the sampling radius increases, together with **second-order descriptors** (two-point covariance and the corresponding radial spectrum) computed within an inscribed circular region.

---

## MATLAB scripts: algorithms, workflow, and dependencies

### Function-designation convention

Each function call is labeled as one of:

- **native:** available in base MATLAB (no specialized toolbox required)
- **native_IPT:** available in MATLAB’s Image Processing Toolbox (IPT)
- **project:** custom function distributed with this project

---

## `playground_detrending_field.m`

### Purpose

Perform an axial stationarity-oriented diagnostic on a 3D cylindrical volume by:

1. Binarizing and masking the volume
2. Computing a slice-wise phase fraction signal **φ(k)** along the **z** direction
3. Detrending **φ(k)** using moving-average windows of varying size
4. Evaluating the detrended residual via its **excess kurtosis**

### Algorithm: axial detrending and excess-kurtosis sweep

**Inputs**
- 3D TIFF volume file: `rockname.tif`
- Window range: 𝓦 = {3, …, 100}

**Outputs**
- Slice-wise phase fraction signal: **φ(k)**
- Excess kurtosis curve: **K_ex(w)**
- Exported figures: `kurtosis.pdf`, `MA.pdf`

**Procedure**

1. **Read volume:**  
   `A <- tiffreadVolume(...)` *(native_IPT)*

2. **Normalize intensities:**  
   `A_d <- mat2gray(A)` *(native_IPT)*

3. **Compute global threshold:**  
   `T <- graythresh(A_d(:))` *(native_IPT)*

4. **Binarize:**  
   `B <- (A_d > T)` *(native)*

5. **Mask cylinder:**  
   `B <- cylinder_masked(B)` *(project)*

6. **Store voxel geometry:**  
   `(N_x, N_y, N_z) <- size(B)` *(native)*

7. **Compute slice-wise phase fraction:** for `k = 1 … N_z`  
   - `x(k) <- k` *(native)*  
   - `φ(k) <- mean(B(:,:,k)(:))` *(native)*

8. **Sweep rolling-window detrending:** for each `w ∈ 𝓦`  
   - `K_ex(w) <- detrend_window_xy(x, φ, w)` *(project)*

9. **Export diagnostic curve:**  
   Plot `K_ex(w)` vs `w` and save as `kurtosis.pdf` *(native)*

10. **Select window:**  
   Choose `w_opt` *(criterion-dependent; native)*

11. **Export detrending visualization:**  
   Run `detrend_window_xy(x, φ, w_opt)` and save as `MA.pdf` *(project + native)*

### Graphical outputs

- `kurtosis.pdf`: excess kurtosis vs moving-average window size *(w vs K_ex)*
- `MA.pdf`: detrending visualization produced by `detrend_window_xy` at the selected window `w_opt`

### Required project functions

- `cylinder_masked.m`
- `detrend_window_xy.m`

---

## `playground_optimal_REV_radius.m`

### Purpose

Compute slice-wise spatial diagnostics from a single 2D cross-sectional slice extracted from a 3D cylindrical volume. The script computes:

1. Disk-averaged phase fraction **⟨φ⟩(r)**
2. Covariance **C(r)** within an inscribed circular region
3. Corresponding radial spectrum **Ĉ(k_r)**
4. A single **2×2** figure summarizing the results

### Algorithm: slice-wise radial averaging, covariance, and spectral diagnostics

**Inputs**
- 3D TIFF file: `filename`
- Slice index: `z0`
- Radius list: `rList`
- Parameters: `doCircleCrop`, `fracR`, `lagStep`, `Nk`

**Outputs**
- Disk-average curve: **⟨φ⟩(r)**
- Covariance: **C(r)**
- Spectrum: **Ĉ(k_r)**
- Exported figure: `REV_r.pdf`

**Procedure**

1. **Read volume:**  
   `V <- tiffreadVolume(filename)` *(native_IPT)*

2. **Binarize and mask:**  
   - `V_n <- mat2gray(V)` *(native_IPT)*  
   - `B <- imbinarize(V_n)` *(native_IPT)*  
   - `B <- cylinder_masked(B)` *(project)*

3. **Extract slice and set geometry:**  
   - `B_slice <- B(:,:,z0)` *(native)*  
   - `(N_x, N_y) <- size(B_slice)` *(native)*  
   - `R_max <- floor(min(N_x, N_y)/2)` *(native)*

4. **Compute disk-averaged phase fraction:** for each `r ∈ rList`  
   - `⟨φ⟩(r) <- circleMeanPlot(B_slice, r, false)` *(project)*

5. **Compute 2D covariance:**  
   - `(rvals, C_r) <- slice2cov_bool(B_slice, doCircleCrop, fracR, lagStep)` *(project)*

6. **Compute radial spectrum:**  
   - `(k_r, Ĉ) <- radialSpectrum2D(rvals, C_r, Nk)` *(project)*

7. **Assemble 2×2 tiled figure:**  
   Use `tiledlayout` / `nexttile` to plot:  
   - (a) slice with overlay circles  
   - (b) `⟨φ⟩(r)` vs radius  
   - (c) `C(r)` vs lag distance  
   - (d) `Ĉ(k_r)` vs `k_r` (log–log)  
   *(native)*

8. **Export:**  
   `print(fig,'REV_r.pdf','-dpdf','-painters')` *(native)*

### Graphical output

- `REV_r.pdf`: **2×2** tiled figure containing  
  (a) binarized slice with overlay circles at `r = r_overlay` and `r = R_max`  
  (b) `⟨φ⟩(r)` vs radius  
  (c) covariance `C(r)` vs lag distance  
  (d) radial spectrum `Ĉ(k_r)` vs `k_r` (log–log)

### Required project functions

- `cylinder_masked.m`
- `circleMeanPlot.m`
- `slice2cov_bool.m`
- `radialSpectrum2D.m`

---

## Dependencies

- **MATLAB** (base)
- **Image Processing Toolbox (IPT)**  
  Required for:
  - `tiffreadVolume`
  - `mat2gray`
  - `graythresh`
  - `imbinarize`

---

## Expected repository layout (suggested)

```text
.
├── playground_detrending_field.m
├── playground_optimal_REV_radius.m
├── cylinder_masked.m
├── detrend_window_xy.m
├── circleMeanPlot.m
├── slice2cov_bool.m
├── radialSpectrum2D.m
└── data/
    └── example_volume.tif
