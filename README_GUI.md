# ULM Master GUI v2.0
### An Interactive Optimization Platform for Ultrasound Localization Microscopy

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-blue)](https://www.mathworks.com/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

---

## 👤 Authors & Affiliation

**Grigori Shapiro**  
Tali Ilovitsh Lab  
School of Biomedical Engineering  
Tel Aviv University, Tel Aviv, Israel

**Supervisor:** Prof. Tali Ilovitsh

---

## 📄 Citation

If you use this software in your research, please cite:

**Thesis:**
> Shapiro, G. "Ultrasound Localization Microscopy with Micro- and Nanobubbles: Processing Framework and Experimental Validation." M.Sc. Thesis, School of Biomedical Engineering, Tel Aviv University, March 2026.

**Journal Article** *(in preparation — citation will be updated upon publication):*
> Shapiro, G., et al. "[Title TBD]." *[Journal TBD]*, [Year TBD].  
> DOI: *pending*

---

## 📋 Table of Contents

1. [Overview](#1-overview)
2. [Scientific Background](#2-scientific-background)
3. [System Requirements](#3-system-requirements)
4. [Installation](#4-installation)
5. [File Structure](#5-file-structure)
6. [Architecture](#6-architecture)
7. [Launching the GUI](#7-launching-the-gui)
8. [Input Data Format](#8-input-data-format)
9. [GUI Layout](#9-gui-layout)
10. [Workflow — Tab by Tab](#10-workflow--tab-by-tab)
    - [Tab 1: Filter](#tab-1-filter)
    - [Tab 2: Detect](#tab-2-detect)
    - [Tab 3: Localize](#tab-3-localize)
    - [Tab 4: Track](#tab-4-track)
    - [Tab 5: Post-Process](#tab-5-post-process)
    - [Tab 6: Render](#tab-6-render)
11. [Parameter Reference](#11-parameter-reference)
12. [Session Management](#12-session-management)
13. [Undo / Redo System](#13-undo--redo-system)
14. [ROI & Masking Tools](#14-roi--masking-tools)
15. [Keyboard Shortcuts](#15-keyboard-shortcuts)
16. [Performance Guide](#16-performance-guide)
17. [Troubleshooting](#17-troubleshooting)
18. [License](#18-license)

---

## 1. Overview

The **ULM Master GUI** is a MATLAB-based interactive platform for end-to-end Ultrasound Localization Microscopy (ULM) processing. It was developed to address a fundamental challenge in ULM research: the high sensitivity required to detect and track sub-micron nanobubbles (NB), which have a substantially weaker acoustic backscatter signal compared to standard microbubbles (MB).

Unlike static "black-box" batch scripts where all parameters are set once and applied globally, the ULM Master GUI provides a **real-time feedback loop**: every parameter adjustment triggers an immediate re-calculation on the cached dataset, visualizing the result instantly. This interactive approach is essential for the systematic optimization of the full processing chain — from clutter filter thresholds to sub-pixel localization constraints to tracking gate sizes — and is particularly critical when working with agents operating near the detection limit.

**Key design principles:**
- Every processing stage is independently tunable with live visual feedback
- The cached SVD decomposition allows instant slider adjustments without re-computation
- Optimized parameters can be directly exported to the batch processing pipeline (`run_ULM_Analysis_Kidney.m`)
- Full session save/load ensures reproducibility across experiments

---

## 2. Scientific Background

ULM is a super-resolution imaging technique that reconstructs microvascular architecture by localizing individual gas-filled contrast agents with sub-pixel precision and linking their positions into trajectories over thousands of frames. Resolution is not limited by the acoustic diffraction limit (typically ~λ/2, hundreds of μm for clinical transducers) but by localization precision, which follows:

```
σ_loc ∝ FWHM / √SNR
```

The pipeline consists of five sequential stages:

| Stage | Challenge | GUI Module |
|-------|-----------|------------|
| **Clutter Filtering** | Separating bubble signal from tissue background | Tab 1: Filter |
| **Detection** | Identifying candidate bubble peaks above noise | Tab 2: Detect |
| **Localization** | Estimating sub-pixel bubble centroid | Tab 3: Localize |
| **Tracking** | Linking localizations into coherent trajectories | Tab 4: Track |
| **Rendering** | Projecting trajectories onto high-resolution grid | Tabs 5–6 |

This GUI was validated on in vitro gelatin phantoms (100–500 μm channels), in vivo rat brain data (benchmark dataset, Chavignon et al.), and in vivo rat kidney data. Full experimental details are available in the associated thesis (see Citation).

---

## 3. System Requirements

### MATLAB Version
- **Minimum:** MATLAB R2020b
- **Recommended:** MATLAB R2022a or later

### Required Toolboxes
| Toolbox | Used For |
|---------|----------|
| Image Processing Toolbox | Masking, ROI tools, morphological operations |
| Signal Processing Toolbox | Butterworth filter, Savitzky-Golay smoothing |
| Statistics and Machine Learning Toolbox | K-means clustering in DCC-SVD |
| Optimization Toolbox | Gaussian fitting (lsqcurvefit) |

### Optional Toolboxes
| Toolbox | Benefit |
|---------|---------|
| Parallel Computing Toolbox | Significant speedup in Gaussian fitting and post-processing via `parfor` |

### Hardware
- **RAM:** 16 GB minimum; 32 GB+ recommended for in vivo datasets
- **Display:** 1920×1080 minimum resolution (GUI is designed for 1600×1000 px)

### Operating System
- Windows 10/11 (primary development platform)
- macOS (tested)
- Linux (tested)

---

## 4. Installation

1. **Download** or clone all files into a single directory:
   ```bash
   git clone https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
   cd YOUR_REPO_NAME
   ```

2. **Open MATLAB** and navigate to the cloned directory.

3. **Add to path** — the GUI does this automatically on launch:
   ```matlab
   addpath(genpath('path/to/ulm-codebase'));
   ```
   Or let the GUI handle it by running from inside the folder.

4. **Verify toolboxes** are installed by running:
   ```matlab
   ver
   ```
   Check that the toolboxes listed in Section 3 are present.

---

## 5. File Structure

```
ulm-codebase/
│
├── [GUI/]
│   ├── ULM_Master_GUI_v2.m       ← Main GUI application (entry point)
│   ├── ULM_Constants.m           ← All default values and limits (centralized)
│   ├── SessionManager.m          ← Save/load session logic
│   ├── UndoRedoManager.m         ← Parameter history management
│   ├── DisplayManager.m          ← All visualization and rendering code
│   └── DataHash.m                ← MD5-based change detection for SVD caching
|
├── [core/]
|   |
│   ├── ULM_Processor.m       ← Main ULM class
|   |
|   ├── [detection/]
|   │   ├── detectBubbles.m       ← Intensity-based regional maxima detection
|   │   ├── detectBubbles_NCC.m   ← Normalized cross-correlation detection
|   │   └── detectBubbles_NP.m    ← Neyman-Pearson hypothesis test detection
|   │
|   ├── [localization/]
|   │   ├── localizeRadialSymmetry.m  ← Gradient-based radial symmetry
|   │   ├── fit2DGaussian.m       ← Full 2D Gaussian NLLS fit
|   │   └── fit2DGaussian_Fast.m  ← Vectorized fast Gaussian fitting
|   │
|   ├── [tracking/]
|   │   ├── trackNearestNeighbor.m    ← Greedy nearest-neighbor linker
|   │   ├── trackHungarian.m          ← Global Hungarian assignment
|   │   ├── trackKalman.m             ← Standard Kalman filter tracker
|   │   ├── trackKalman_Advanced.m    ← Hierarchical Kalman tracker (HKT)
|   │   ├── calculateCostMatrix.m     ← Smart Cost Matrix (SCM) engine
|   │   └── munkres.m                 ← Munkres/Hungarian solver
|   │
|   ├── [filtering/]
|   │   ├── SVD_filter.m          ← Standard SVD clutter filter
|   │   ├── SVD_SSM.m             ← Spatial Similarity Matrix thresholding
|   │   ├── SVD_blockwise.m       ← Block-wise adaptive SVD
|   │   ├── DCC_SVD.m             ← Density Canopy Clustering SVD
|   │   ├── run_SVD_Decomposition.m
|   │   ├── reconstruct_SVD_Signal.m
|   │   └── Butterworth_bandpass_filter.m
|   │
|   └── [rendering/]
|       ├── renderHistogram.m     ← Histogram accumulation (default)
|       └── renderGaussian.m      ← Gaussian-blurred density map
|
├── [config/]
│   ├── setDefaultParams.m        ← Experiment parameter configuration
│   └── getExpParams.m            ← Parses info.txt metadata file
│
├── [utils/]
│   ├── applyQualityControl.m
│   ├── applyAccelerationConstraint.m
│   ├── applyDirectionConstraint.m
│   ├── applyVDConstraint.m
│   ├── generateVesselMask.m
│   ├── printTrackMetrics.m
│   └── analyze_ULM_Features.m
│
└── [examples/]
    ├── run_ULM_Analysis_Kidney.m  ← Full batch pipeline script
    ├── split_ImageData_tot_Kidney.m
    ├── Bmode_video.m
    └── Kidney_video.m
```

---

## 6. Architecture

The GUI follows a strict **separation of concerns** with four manager classes and one monolithic main file:

```
ULM_Master_GUI_v2.m          (Layout builder + all callbacks)
       │
       ├── ULM_Constants.m   (Read-only: all defaults and limits)
       ├── SessionManager.m  (Save/load: serializes the full app struct)
       ├── UndoRedoManager.m (History: up to 20 parameter states)
       └── DisplayManager.m  (Visualization: frame rendering, overlays)
```

### State Machine
The GUI state is tracked in `app.state.currentState`, which controls which panels are enabled. States advance sequentially:

```
-1 → No data loaded    (only Load button active)
 0 → Data loaded       (Filter tab unlocked)
 1 → Filtered          (Detect tab unlocked)
 2 → Detected          (Localize tab unlocked)
 3 → Localized         (Track tab unlocked)
 4 → Tracked           (Post-Process tab unlocked)
 5 → Post-processed    (Render tab unlocked)
```

### SVD Caching
SVD decomposition (the most expensive step) is cached using a data hash (`DataHash.m`). If the raw data and cutoff range have not changed since the last computation, the cached `U`, `S`, `V` matrices are reused. This makes slider adjustments for SVD cutoff **instantaneous**.

### Debounced Display
Display updates triggered by sliders are debounced with a 50 ms delay (configurable in `ULM_Constants.DEBOUNCE_DELAY`) to prevent UI freezing during rapid slider movement.

---

## 7. Launching the GUI

```matlab
% Simply run from the MATLAB command window:
ULM_Master_GUI_v2
```

The GUI will:
1. Attempt to load default parameters via `setDefaultParams.m`
2. Fall back to safe system defaults if `setDefaultParams.m` is missing or fails
3. Display a warning dialog if the fallback was used
4. Open a 1600×1000 px figure window

---

## 8. Input Data Format

### Required Format
- **File type:** `.mat` (MATLAB data file)
- **Variable:** Any 3D numeric array of shape `[Nz × Nx × Nt]`
  - `Nz` = number of axial pixels (image height)
  - `Nx` = number of lateral pixels (image width)
  - `Nt` = number of frames (time dimension)
- The GUI automatically detects 3D numeric arrays in the loaded `.mat` file regardless of the variable name.

### Physical Calibration
Physical pixel sizes are read from `info.txt` (via `getExpParams.m`) or set manually in the GUI's **Fundamental Parameters** section:

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `Framerate (Hz)` | Acquisition frame rate | 200–1000 Hz |
| `Pixel Size X (mm)` | Lateral pixel size | ~0.05–0.1 mm |
| `Pixel Size Z (mm)` | Axial pixel size | ~0.05–0.1 mm |

These values affect all downstream physical calculations (velocity in mm/s, linking distances in mm).

### Example: Loading Data
```matlab
% Data must be a 3D array in a .mat file:
data = rand(122, 260, 1500);   % [Nz x Nx x Nt]
save('my_data.mat', 'data');
% Then click "Load Raw Data (.mat)" in the GUI
```

---

## 9. GUI Layout

The interface is divided into three functional zones:

```
┌─────────────────────────────────────────────────────────────────┐
│  Menu Bar: [Load Session] [Save Session] [Undo] [Redo]  [Memory]│
├────────────────────────────────────────┬────────────────────────┤
│                                        │  Workflow Tabs:        │
│                                        │  1.Filter 2.Detect ... │
│                                        │  ──────────────────    │
│   Visualization Canvas                 │  Parameter             │
│   (Central, real-time display)         │  Control Panel         │
│                                        │  (Active tab controls) │
│                                        │                        │
│                                        │  [Run Button]          │
│                                        │  [Reset to Defaults]   │
├────────────────────────────────────────┴────────────────────────┤
│  Frame Slider ────────────────────────────────── [Frame: ___]   │
│  [Load Raw Data]           [▶ Play / ■ Stop]       Status ●     │
└─────────────────────────────────────────────────────────────────┘
```

### Menu Bar
| Button | Function |
|--------|----------|
| 📁 Load Work Session | Restores a previously saved complete GUI state |
| 💾 Save Work Session | Saves all data, parameters, and results to a `.mat` file |
| ↶ Undo | Reverts last parameter change (up to 20 levels) |
| ↷ Redo | Re-applies an undone change |
| Memory: X MB | Real-time RAM usage monitor (auto-updates every 2 s) |

### Status Bar
The status lamp and label in the bottom-right indicate the current state:
- 🟢 **Green** — Ready / operation complete
- 🔴 **Red** — Processing in progress
- 🔵 **Blue** — Informational (e.g., "Undo complete")

---

## 10. Workflow — Tab by Tab

### Tab 1: Filter

**Purpose:** Remove stationary tissue clutter (high-energy, temporally coherent) from the IQ data, retaining the dynamic microbubble signal.

#### Step A: Spatial Crop (Optional)
Before filtering, you can spatially crop the dataset to a sub-region of interest. This significantly reduces SVD computation time.

| Parameter | Description |
|-----------|-------------|
| **Enable Spatial Crop** | Master on/off switch |
| **Draw New Crop** | Interactively draw a crop rectangle on the visualization canvas |
| **Load Existing Crop** | Load a previously saved `cropBox.mat` |
| **Save Current Crop** | Export the current crop rectangle for reuse |

> **Why crop?** For a 122×260 frame, SVD operates on a 31,720×Nt matrix. Cropping to a 60×130 ROI reduces this by 4×, giving a ~16× speedup.

#### Step B: Clutter Filter Method

Select one of four available methods from the **Filter Method** dropdown:

---

**`svd_filter` — Standard SVD Filter**

The baseline method. Decomposes the Casorati matrix `X = UΣV*` and retains only singular values within the specified range `[Cutoff Start, Cutoff End]`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Cutoff Start** | 4 | First singular value index to keep. Values below this (indices 1 to Cutoff Start-1) are discarded as tissue clutter. Increase to remove more tissue. |
| **Cutoff End** | 450 | Last singular value index to keep. Values above this are discarded as noise. Decrease to remove noise; increase to retain more signal. |

> **Caching:** The full SVD decomposition (`U`, `Σ`, `V`) is computed once and cached. Subsequent slider adjustments reconstruct the filtered data instantly from the cache without recomputation.

> **Tuning guide:** Start with a broad range (e.g., [4, 1000]). Gradually increase Cutoff Start until tissue clutter disappears. Gradually decrease Cutoff End until background noise is suppressed. Optimal range isolates moving bubbles.

---

**`svd_ssm` — SVD with Spatial Similarity Matrix (SSM)**

Implements the method of Baranger et al. [2023]. Calculates the Pearson correlation between absolute values of spatial singular vectors (`U`) to identify the transition boundary between tissue-dominated and bubble-dominated components automatically.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Cutoff Start** | 4 | Lower bound (same interpretation as standard SVD) |
| **Cutoff End** | 450 | Upper bound; the SSM provides a statistical guide for this |

> **When to use:** When manual SVD thresholding is ambiguous. The SSM correlation matrix provides an objective visual cutoff.

---

**`dcc_svd` — Density Canopy Clustering SVD**

Implements the DCC-SVD method of Han et al. [2024]. Each singular component is characterized by a 3D feature vector:
1. Normalized log-energy
2. Power-weighted central temporal frequency
3. Spatial correlation to the mean spatial vector

K-means clustering (seeded by density canopy centers) then partitions components into "Tissue", "Blood", and "Noise" clusters.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Cutoff Start** | 1 | Lower bound before DCC classification |
| **DCC Components** | Auto | Number of cluster representatives shown |

> **When to use:** For fully automated thresholding without manual SVD tuning. Particularly effective when tissue and bubble signatures are spectrally similar.

---

**`svd_blockwise` — Block-Wise Adaptive SVD**

Divides the image into overlapping spatial blocks and applies independent adaptive SVD thresholding to each block. Accounts for spatially varying clutter characteristics.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Block Size (mm)** | 4.0 mm | Spatial extent of each processing block (square). Set to `[]` for automatic selection. |
| **Overlap (%)** | 75% | Block overlap percentage. Higher overlap = better reconstruction quality but slower computation. Use 75% for exploratory runs, 93.75% for publication-quality results. |
| **Threshold Method** | DopplerGradient | Strategy for determining per-block SVD cutoff: `DopplerGradient`, `SSM`, `Hybrid`, or `Manual` |
| **Tissue Freq (Hz)** | -1 (auto) | Tissue Doppler frequency threshold. Set to -1 for automatic estimation: `max(5, min(20, framerate/50))`. Increase (e.g., 15 Hz) for fast-moving tissue. |
| **MP Deviation σ** | 2.0 | Marchenko-Pastur sensitivity for the high cutoff. Higher value = fewer components classified as blood. |
| **Gradient Pct** | 0.10 | Sensitivity of inflection detection for Cutoff 1A. Lower = earlier inflection detection. |
| **Min Blood Comps** | 3 | Minimum number of blood components per block (floor constraint). |
| **Max Tissue Frac** | 0.60 | Maximum fraction of singular values classifiable as tissue (ceiling constraint). |
| **Plot Maps** | false | If true, displays spatial threshold maps after filtering. |

---

#### Step C: Butterworth Bandpass Filter (Optional)

An auxiliary temporal frequency filter applied after SVD to further isolate the bubble signal within a specific frequency band.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Enable** | Off | Master switch |
| **Cutoff [Low, High] (Hz)** | [50, 250] | Passband frequencies. Tune based on expected bubble velocity and acquisition framerate. |
| **Order** | 2 | Filter order. Higher order = sharper rolloff. |

> **When to use:** Primarily for high-frame-rate datasets (≥1000 Hz) where specific frequency isolation is needed, or as a supplement to SVD when slow-moving bubble signals partially overlap with tissue.

#### Step D: Spatial Filter (Optional)

A per-frame spatial convolution filter for additional noise reduction.

| Parameter | Options | Description |
|-----------|---------|-------------|
| **Method** | `None`, `Gaussian`, `Median`, `DoG`, `Top-Hat` | Type of spatial filter |
| **Kernel Size (px)** | Default: 3 | Filter kernel size (must be odd) |
| **Sigma 1** | Default: 0.5 | Primary Gaussian sigma |
| **Sigma 2** | Default: 4.0 | Secondary sigma (for Difference-of-Gaussians) |

#### Running the Filter
Click **"Run Filter"**. The status lamp turns red during processing and green upon completion. The visualization canvas updates to show the filtered frame.

---

### Tab 2: Detect

**Purpose:** Identify candidate bubble locations in the filtered frames as integer-pixel regional maxima.

#### Step A: ROI Masking (Vessel Map)

Restricts detection to anatomically relevant regions, reducing false positives and computation.

**Manual ROI:**

| Control | Description |
|---------|-------------|
| **Load Mask** | Load a pre-existing binary mask (`Mask.mat`) |
| **Create New Mask** | Interactively draw a polygon mask on the canvas |
| **Save Current Mask** | Export the current mask to `Mask.mat` |
| **Status: None/Loaded** | Indicator showing current mask state |

**Algorithmic Vessel Masking (Auto-Mask):**

Generates a "Flow Probability Mask" from the temporal average of the SVD-filtered sequence.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Enable** | Off | Master switch for algorithmic masking |
| **Enhancement Method** | `Top-Hat` | `None`, `CLAHE`, `Top-Hat` (vesselness), or `Sharpen` |
| **Strength** | 1.0 | Enhancement intensity (0 to 1) |
| **Gamma** | 1.0 | Contrast adjustment (< 1 = brighter, > 1 = darker) |
| **Threshold** | 0.081 | Binary threshold applied to the enhanced temporal average image |

**Enhancement methods explained:**

| Method | Best For | Effect of Strength Slider |
|--------|----------|--------------------------|
| **CLAHE** | Uniform regions with low local contrast | Controls clip limit (higher = more contrast enhancement) |
| **Top-Hat** | Vessel-like elongated structures | Controls structural element radius (vessel width estimate) |
| **Sharpen** | Improving edge definition | Controls sharpening intensity |

#### Step B: Detection Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| **Detection Method** | `Intensity` | `Intensity`, `NP`, `NCC` | Algorithm for finding candidate peaks |
| **Intensity Threshold** | 0.15 | 0.01–1.0 | Normalized intensity threshold. Candidates below this fraction of the frame's maximum are rejected. Lower = more candidates (higher sensitivity), higher = fewer candidates (higher specificity). |
| **Max Bubbles/Frame** | 100 | 1–5000 | Maximum number of candidates retained per frame. **Critical:** if the detected count equals this limit, the algorithm is saturated and valid bubbles are being missed. Increase until the detected count drops below the limit. |

**Detection Methods:**

| Method | Description | Key Parameter |
|--------|-------------|---------------|
| **Intensity** | Regional maxima above a normalized intensity threshold. Fast and robust. | `Intensity Threshold` |
| **NP (Neyman-Pearson)** | Hypothesis test with controlled false alarm rate. Statistically rigorous. | `NP Alpha` (false alarm rate, e.g., 1e-4) |
| **NCC (Normalized Cross-Correlation)** | Template matching with a reference PSF (Gaussian or experimental). | `NCC Threshold` (minimum correlation, e.g., 0.7), `PSF Type`, `PSF Size` |

**FWHM (Full Width at Half Maximum):**

| Parameter | Default | Description |
|-----------|---------|-------------|
| **FWHM [X, Z] (px)** | [3, 3] | Estimated PSF size in pixels. Defines the ROI size around each candidate for the subsequent localization step. Should approximate the bubble's appearance in the filtered image. |

#### Running Detection
Click **"Run Detection (Masked)"**. Detected candidates are overlaid as red markers on the visualization canvas. The status bar displays the count per frame (e.g., "Detection: 141 Particles | Frame 1/1500").

> **Saturation check:** If the status reads exactly `Max Bubbles/Frame` candidates, increase the limit.

---

### Tab 3: Localize

**Purpose:** Refine integer-pixel candidate positions to sub-pixel precision using the PSF geometry.

#### Localization Method

| Method | Algorithm | Speed | Precision | Best For |
|--------|-----------|-------|-----------|---------|
| **`radial`** | Gradient-based radial symmetry (Parthasarathy, 2012) | Fast (non-iterative) | High | High-throughput, moderate SNR |
| **`gaussian_fit`** | 2D Gaussian NLLS (Levenberg-Marquardt) | Slow (iterative) | Highest (Gold Standard) | Low-count, high-SNR data |
| **`gaussian_fit_fast`** | Vectorized Gaussian fitting with parallel processing | Medium | High | Large datasets with Parallel Computing Toolbox |

The sub-pixel center `(x̂, ŷ)` for Radial Symmetry is the weighted least-squares intersection of gradient lines:

```
weights = |∇I|²   (pixels on PSF slopes carry most information)
```

The Gaussian model is:
```
I(x,y) = A · exp(−((x−x₀)² + (y−y₀)²) / 2σ²) + C
```

#### Localization Quality Control (QC)

A multi-layer filter removes candidates that do not conform to the expected PSF shape:

| QC Check | Default | Description |
|----------|---------|-------------|
| **Divergence Check** | On | Rejects sub-pixel solutions that shift more than `Max Shift Factor × FWHM/2` pixels from the coarse integer peak. Prevents convergence failures. |
| **Max Shift Factor** | 1.0 | Multiplier on FWHM/2. Increase to allow larger shifts (for noisy data). |
| **ROI Maxima Check** | On | Rejects ROIs with multiple intensity peaks (overlapping bubbles). |
| **Max ROI Maxima** | 3 | Maximum number of local maxima allowed in an ROI before rejection. |
| **Min Gradient²** | 1e-6 | Minimum gradient magnitude for radial symmetry (rejects flat regions). |
| **Min Determinant** | 1e-6 | Minimum matrix determinant for radial symmetry linear system (rejects numerically unstable fits). |
| **R² Threshold** | 0.3 | Minimum goodness-of-fit for Gaussian fitting (rejects poor fits). |
| **Gaussian Box Radius** | 2 px | Half-size of the ROI used for Gaussian fitting. |

#### Running Localization
Click **"Run Localization"**. The status bar shows the localization yield (e.g., "Localization: 107 Particles | Frame 1/1500"), i.e., how many of the detected candidates survived QC.

> **Yield monitoring:** A drop from 141 detections to 107 localizations (24% rejection) is typical and healthy. A rejection rate >60% suggests the FWHM or QC parameters need adjustment.

---

### Tab 4: Track

**Purpose:** Link sub-pixel localizations across frames into continuous trajectories.

#### Tracking Algorithm

| Algorithm | File | Description |
|-----------|------|-------------|
| **`nn` (Nearest Neighbor)** | `trackNearestNeighbor.m` | Greedy local linker. Fast but prone to fragmentation in high-density regions. |
| **`Hungarian` (HT)** | `trackHungarian.m` | Global linear assignment (Munkres solver). Improves directionality vs. NN but still lacks predictive capability. |
| **`Kalman` (KT)** | `trackKalman.m` | Predictive Kalman filter with constant-velocity motion model. Bridges temporal gaps ("blinking"). Recommended for most datasets. |
| **`Kalman_Advanced` (HKT)** | `trackKalman_Advanced.m` | Hierarchical multi-pass tracker (Taghavi et al., 2022). Processes velocity ranges sequentially (slow → fast). Best for heterogeneous flow (e.g., kidney: arcuate arteries + peritubular capillaries simultaneously). |

#### Core Tracking Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| **Max Linking Distance (px)** | 2 | 0.1–10 | Maximum Euclidean distance allowed for linking a localization to an active track. In Kalman mode, this is the distance from the *predicted* state. |
| **Max Gap Closing (frames)** | 1 | 0–10 | Number of consecutive missing frames a track can survive before being terminated. Set to 0 to disable gap closing. |
| **Min Track Length** | 5 | 2–20 | Minimum number of localizations for a track to be retained. Shorter tracks are discarded as fragments or noise. |

#### Kalman Filter Settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Motion Model** | `ConstantVelocity` | State vector `[x, y, vx, vy]`. Use `ConstantAcceleration` for rapidly changing velocity. |
| **Process Noise** | 3 | Kalman model flexibility. Low = tracker trusts the motion model (stiff). High = tracker follows measurements more closely (flexible). Increase if tracks break on sharp turns. |
| **Assignment Method** | `hungarian` | Inner assignment solver: `hungarian` (global, optimal) or `nn` (nearest neighbor, fast). |

#### Smart Cost Matrix (SCM)

When **Use Advanced Cost Matrix** is enabled, the linking cost is computed as:

```
C_total(i,j) = C_dist · (1 + W_slope · P_angle) · (1 + W_int · P_intensity)
```

Where:
- `C_dist` = base spatial distance (Euclidean or Kalman-predicted)
- `P_angle = max(0, θᵢⱼ − θ_gate)` = directional penalty for turns beyond the safety gate
- `P_intensity = |I_current − Ī_track| / (Ī_track + ε)` = brightness consistency penalty

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Direction Weight (W_dir)** | 2 | Penalty weight for direction changes (0 = disabled, typical 1–5) |
| **Angle Penalty Slope (W_slope)** | 0.3 | Linear slope of the directional cost above the gate angle |
| **Brightness Weight (W_int)** | 2 | Penalty weight for brightness changes (0 = disabled, typical 0–3) |
| **Max Angle Gate (°)** | 70 | Soft limit: angles above this incur an increasing penalty |
| **Hard Angle Gate (°)** | 90 | Hard limit: angles above this set cost to infinity (physically impossible links are blocked) |
| **Direction History** | 4 pts | Number of past positions used to estimate current trajectory direction |

#### Hierarchical Kalman Tracker (HKT) Settings

The HKT decomposes tracking into `N` velocity levels, processing slow bubbles first and subtracting their localizations before tracking faster ones.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Max Velocity (mm/s)** | 20 | Global upper velocity limit across all levels |
| **Num Levels** | 5 | Number of velocity bands (e.g., 5 levels over 0–20 mm/s) |
| **Spacing Power** | 1.0 | 1.0 = linear spacing. >1.0 = more levels at low velocities. <1.0 = more at high velocities. |
| **Enable Overlap** | On | Add an overlap band between adjacent velocity levels to prevent missed assignments at boundaries |
| **Overlap Width (mm/s)** | 2.0 | Width of the overlap band |
| **Process Noise α** | 0.01 | Process noise multiplier: `σ_process = α × v_max_level`. Increase if tracks break on turns. |
| **Measurement Noise β** | 0.025 | Measurement noise base: `σ_meas = β / 2^(level−1)`. Decrease if localizations are highly precise. |
| **Forward-Backward** | On | Dual-pass tracking: runs forward then backward through time, maximizing track yield |

#### Post-Tracking Quality Control

Optional filters applied to *completed* trajectories:

| QC Constraint | Default | Description |
|--------------|---------|-------------|
| **Direction Constraint** | Off | Rejects tracks with turns larger than `Max Angle (°)`. Use for phantoms with known straight channels. |
| **Acceleration Constraint** | Off | Adaptive acceleration gating. Rejects unrealistic velocity jumps. |
| **Velocity Dispersion (VD) Constraint** | Off | Rejects "jittery" tracks where path length >> displacement (high tortuosity). |
| **Max VD Ratio** | 3.0 | Maximum allowed ratio of path length to net displacement (Tortuosity Index threshold). |

#### Running Tracking
Click **"Run Tracking"**. Tracks are overlaid on a Temporal Mean Intensity Projection (TMIP) of the filtered data. The track count and mean tortuosity are printed to the MATLAB console.

---

### Tab 5: Post-Process

**Purpose:** Smooth discrete trajectory points and interpolate onto a sub-pixel grid before rendering.

#### Track Smoothing

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| **Enable Smoothing** | On | On/Off | Master switch |
| **Window Size** | 11 | 3–21 (odd) | Smoothing window width. Larger = smoother but may lose fine vessel curvature. |
| **Smoothing Method** | `sgolay` | `sgolay`, `rloess`, `gaussian`, `movmean` | Algorithm for smoothing raw track points. Savitzky-Golay (`sgolay`) is recommended: it suppresses high-frequency jitter while preserving peak positions and vessel geometry. |

**Smoothing method comparison:**

| Method | Preserves Shape | Noise Suppression | Speed |
|--------|----------------|------------------|-------|
| `sgolay` (Savitzky-Golay) | ★★★★★ | ★★★★ | Fast |
| `gaussian` | ★★★★ | ★★★★★ | Fast |
| `movmean` (Moving Average) | ★★★ | ★★★ | Fastest |
| `rloess` (Robust Loess) | ★★★★ | ★★★ | Slow |

#### Track Interpolation

After smoothing, tracks are densely interpolated to ensure continuity on the super-resolution grid.

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| **Interpolation Method** | `spline` | `spline`, `pchip`, `linear`, `makima` | Interpolation algorithm. `spline` is recommended for smooth vessel profiles. |
| **Interpolation Step** | 1/upsampling_factor | — | Sub-pixel spacing between interpolated points (auto-set from upsampling factor) |

#### Min Track Length Filter (Display)

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Final Min Length** | Adjustable | Slider to set the minimum track length for *display* purposes, without re-running tracking. This allows rapid exploration of the length/density trade-off. |

#### Running Post-Processing
Click **"Run Post-Processing"**. Smoothed tracks are displayed in the canvas.

---

### Tab 6: Render

**Purpose:** Project all trajectories onto a high-resolution grid to generate super-resolution images.

#### Rendering Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| **Upsampling Factor** | 10 | 1–10 | Final image resolution multiplier relative to the native pixel size. Factor 10 on a 50 μm pixel gives 5 μm super-resolution. |
| **Render Method** | `histogram` | `histogram`, `gaussian` | `histogram` = unbiased count accumulation (recommended). `gaussian` = each localization is rendered as a small Gaussian blob (smoother appearance). |
| **Gaussian Sigma** | 0.3 | — | Spread of the Gaussian blob in the `gaussian` render mode (in super-resolved pixels). |

#### Generated Outputs

Clicking **"Generate & Display Final Images"** produces four super-resolution maps simultaneously:

| Map | Description |
|-----|-------------|
| **Density Map** | Accumulates bubble counts per super-resolved pixel. Power-law compression (γ = 0.3) applied to simultaneously visualize high-flux arteries and sparse capillaries. |
| **Raw Velocity Map** | Arithmetic mean of instantaneous velocity vectors (|p_t − p_{t-1}| / Δt) at each pixel. Unbiased statistical representation. |
| **Filtered Velocity Map** | Gaussian-smoothed (σ = 0.6 super-resolved pixels) version of the raw velocity map. Bridges discrete sampling gaps; suppresses isolated high-velocity outliers from tracking errors. |
| **Combined Fusion Map (HSV)** | Dual-mode visualization: **Hue** = velocity (blue → red), **Value** (brightness) = local vessel density. Correlates anatomy with hemodynamics in a single image. |

All maps are saved as `.png` and `.fig` files in the data folder's `Results/` subdirectory.

---

## 11. Parameter Reference

### Default Values Summary (from `ULM_Constants.m`)

| Category | Parameter | Default | Range |
|----------|-----------|---------|-------|
| **Filter** | SVD Cutoff | [1, 100] | [1, Nt] |
| **Filter** | Butterworth Cutoff (Hz) | [10, 100] | (0, framerate/2) |
| **Filter** | Butterworth Order | 4 | 1–8 |
| **Detection** | Intensity Threshold | 0.30 | 0.01–1.0 |
| **Detection** | Max Bubbles/Frame | 2000 | 1–5000 |
| **Detection** | FWHM [X,Z] (px) | [1.5, 1.5] | — |
| **Localization** | Gaussian Box Radius | 3 px | 2–10 |
| **Localization** | Max Shift Factor | 2.0 | — |
| **Tracking** | Max Linking Distance | 2.0 px | 0.1–10 |
| **Tracking** | Max Gap Closing | 2 frames | 0–10 |
| **Tracking** | Min Track Length | 3 | 2–20 |
| **Tracking** | Kalman Process Noise | 10 | — |
| **QC** | Max Angle Change | 90° | — |
| **QC** | Acceleration C Factor | 3.0 | — |
| **QC** | VD Ratio | 2.0 | — |
| **Post-Process** | Smoothing Window | 5 (odd) | 3–21 |
| **Post-Process** | Interpolation Step | 0.5 | — |
| **Render** | Upsampling Factor | 3 | 1–10 |
| **ROI** | Gamma | 1.0 | 0.1–3.0 |
| **ROI** | Threshold | 0.0 | — |
| **UI** | Debounce Delay | 50 ms | — |
| **UI** | Memory Update Interval | 2.0 s | — |
| **Undo** | Max History States | 20 | — |

---

## 12. Session Management

The session system allows complete workspace preservation and reproducibility.

### Saving a Session
Click **💾 Save Work Session** in the menu bar. A `.mat` file is created containing:
- Loaded raw data
- All processed intermediate results (filtered data, detections, localizations, tracks)
- All current parameter values
- Mask and crop settings
- GUI state (current tab, frame index)

### Loading a Session
Click **📁 Load Work Session**. All data and parameters are restored and the GUI updates automatically to the correct processing state.

### Quick Save
Use **Ctrl+S** (if configured) or re-click Save to overwrite the last saved session file without a dialog.

### Session Info (without full load)
From the MATLAB command window:
```matlab
sm = SessionManager();
info = sm.getSessionInfo('my_session.mat');
```
Returns metadata (date, processing state, data size) without loading the full dataset.

---

## 13. Undo / Redo System

The undo system tracks parameter changes with named checkpoints.

| Action | Control | Notes |
|--------|---------|-------|
| **Undo** | ↶ button or **Ctrl+Z** | Reverts last parameter change |
| **Redo** | ↷ button or **Ctrl+Y** | Re-applies undone change |
| **History** | `app.undoManager.displayHistory()` | Prints full undo chain to console |

- Up to **20 levels** of undo are maintained (configurable via `ULM_Constants.MAX_UNDO_STATES`)
- Undo/Redo only tracks **parameter changes**, not processing results
- After undoing, click the relevant "Run" button to re-process with the restored parameters
- The ↶ and ↷ buttons are automatically greyed out when no history is available

---

## 14. ROI & Masking Tools

Masks restrict detection and processing to anatomically meaningful regions, reducing false positives and computation time.

### Manual Polygon Mask
1. In Tab 2, click **"Create New Mask"**
2. Draw a polygon outline on the Visualization Canvas
3. Double-click to close the polygon
4. Click **"Save Current Mask"** to export as `Mask.mat`

### Algorithmic Auto-Mask
The auto-mask generates a binary Flow Probability Mask from the filtered image sequence:
1. Compute temporal mean of SVD-filtered frames
2. Apply selected enhancement (CLAHE / Top-Hat / Sharpen)
3. Apply Gamma contrast adjustment
4. Threshold to produce binary mask

**CLAHE parameters (`ULM_Constants.m`):**
- Clip limit range: 0.001–0.05
- Number of tiles: 8×8

**Top-Hat parameters:**
- Structural element radius range: 1–11 px

### Applying the Mask
With a mask loaded, detection (`detectBubbles.m`) will only evaluate pixels where the mask is non-zero. Pixels outside the mask are ignored.

---

## 15. Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| **Ctrl+Z** | Undo last parameter change |
| **Ctrl+Y** | Redo |
| **Spacebar** | Play / Pause frame playback |
| **← Arrow** | Previous frame (when frame slider is focused) |
| **→ Arrow** | Next frame (when frame slider is focused) |

---

## 16. Performance Guide

### For Large Datasets (>100 MB)

1. **Spatial crop first** — Reduces SVD matrix size quadratically. A 2× crop = 4× SVD speedup.

2. **Use the cached SVD** — After the first filter run, slider adjustments are instantaneous. Avoid clicking "Run Filter" again unless the data or method changes.

3. **Enable Parallel Computing Toolbox** — Speeds up Gaussian fitting and post-processing:
   ```matlab
   parpool('local', 4);   % Start 4 parallel workers
   ```

4. **Limit playback updates** — Disable "Preview ROI Overlay" during parameter tuning.

5. **Use `histogram` render mode** — The `gaussian` render mode is significantly slower for large track sets.

6. **Monitor memory** — The memory monitor in the menu bar updates every 2 seconds. If usage approaches system RAM, clear intermediate variables:
   ```matlab
   app = guidata(gcf);
   app.data.U = [];
   app.data.S_diag = [];
   app.data.V = [];
   guidata(gcf, app);
   ```

### SVD Cutoff Tuning Workflow (Fastest Approach)
1. Run filter once with a wide range (e.g., [1, 1000])
2. Open the frame slider and scrub to a frame with visible bubbles
3. Slowly increase Cutoff Start until tissue background disappears
4. Slowly decrease Cutoff End until noise floor is suppressed
5. All adjustments after the first run are instant (cached SVD)

---

## 17. Troubleshooting

### GUI won't open
- Ensure all `.m` files are in the same directory and on the MATLAB path
- Run from the directory containing `ULM_Master_GUI_v2.m`
- Check that required toolboxes are licensed: `ver`

### "Out of Memory" error during SVD
- Reduce data size by spatial cropping (Tab 1)
- Close other MATLAB figures: `close all`
- Use 64-bit MATLAB (required for >2 GB variables)
- Clear workspace: `clear all`

### Detection count always equals Max Bubbles/Frame
- The detector is saturated. Increase "Max Bubbles/Frame" until the count drops below the limit.

### SVD filtering shows no change when moving sliders
- The SVD cache may be stale. Click "Run Filter" once to recompute and re-cache.
- Verify the displayed frame index is a frame with visible signal.

### Tracking produces too many short fragmented tracks
- Increase "Min Track Length" to filter fragments
- Enable the Smart Cost Matrix (Tab 4)
- Reduce "Max Linking Distance" to prevent cross-linking between nearby vessels
- Switch from NN/Hungarian to Kalman tracking for predictive gap-bridging

### Localization yield is <30% of detections
- FWHM setting may not match the actual bubble PSF size — adjust in Tab 2
- Reduce "Max Shift Factor" if bubbles are drifting too far from peaks
- Check if the filtered image has significant residual noise (re-tune filter cutoffs)

### Display not updating after parameter change
- Check status lamp: if red, processing is still running
- Click the frame number field and press Enter to force a refresh
- If the issue persists, save the session and restart the GUI

### Undo button greyed out
- At least 2 parameter states are needed for undo. Make one parameter change and the button will enable.

### Session fails to load
- Ensure the session file was saved with the same version of the codebase
- Check that all required `.m` files are present on the path

---

## 18. License

This software is released under the **MIT License**.

```
MIT License

Copyright (c) 2026 Grigori Shapiro, Tali Ilovitsh Lab,
Tel Aviv University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
```

---

## Acknowledgments

This work was carried out under the supervision of **Prof. Tali Ilovitsh** at the School of Biomedical Engineering, Tel Aviv University.

The Hierarchical Kalman Tracker (`trackKalman_Advanced.m`) is based on the methodology described in:
> Taghavi, I. et al. "Ultrasound super-resolution imaging with a hierarchical Kalman tracker." *Ultrasonics*, 124, 106742, 2022. DOI: 10.1016/j.ultras.2022.106695

The SSM clutter filtering method (`SVD_SSM.m`) is based on:
> Baranger, J. et al. "Fast Thresholding of SVD Clutter Filter Using the Spatial Similarity Matrix." *IEEE TUFFC*, 70(8), 821–830, 2023. DOI: 10.1109/TUFFC.2023.3289235

The DCC-SVD method (`DCC_SVD.m`) is based on:
> Han, X. et al. "An adaptive spatiotemporal filter for ultrasound localization microscopy based on density canopy clustering." *Ultrasonics*, 144, 107446, 2024. DOI: 10.1016/j.ultras.2024.107446

The benchmark in vivo rat brain dataset used for validation was provided by:
> Chavignon, A. et al. "In vivo rat brain for Ultrasound Localization Microscopy." Zenodo, 2023. DOI: 10.5281/zenodo.7883227

The radial symmetry localization algorithm (`localizeRadialSymmetry.m`) is based on:
> Parthasarathy, R. "Rapid, accurate particle tracking by calculation of radial symmetry centers." *Nature Methods*, 9(7), 724–726, 2012. DOI: 10.1038/nmeth.2071

---

*For questions, bug reports, or feature requests, please open an Issue on the GitHub repository.*
