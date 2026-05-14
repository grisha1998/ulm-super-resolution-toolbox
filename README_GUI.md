# ULM Master GUI v3.0 - Definitive Documentation Manual

### An Interactive Optimization Platform for Ultrasound Localization Microscopy

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-blue)](https://www.mathworks.com/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

---

## Authors & Affiliation

**Grigori Shapiro**
Tali Ilovitsh Lab
School of Biomedical Engineering
Tel Aviv University, Tel Aviv, Israel

**Supervisor:** Prof. Tali Ilovitsh

---

## Citation

If you use this software in your research, please cite:

**Thesis:**
> Shapiro, G. "Ultrasound Localization Microscopy with Micro- and Nanobubbles: Processing Framework and Experimental Validation." M.Sc. Thesis, School of Biomedical Engineering, Tel Aviv University, March 2026.

**Journal Article** *(in preparation — citation will be updated upon publication):*
> Shapiro, G., et al. "[Title TBD]." *[Journal TBD]*, [Year TBD].
> DOI: *pending*

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Scientific Background](#2-scientific-background)
3. [System Requirements](#3-system-requirements)
4. [Installation & Launch](#4-installation--launch)
5. [File Structure](#5-file-structure)
6. [Architectural Deep-Dive](#6-architectural-deep-dive)
    - [6.1 Design Patterns & Separation of Concerns](#61-design-patterns--separation-of-concerns)
    - [6.2 The Central Algorithm Registry](#62-the-central-algorithm-registry)
    - [6.3 Pipeline State Machine & Downstream Invalidation](#63-pipeline-state-machine--downstream-invalidation)
    - [6.4 SVD Caching via DataHash](#64-svd-caching-via-datahash)
    - [6.5 Debounced Display & the DisplayManager](#65-debounced-display--the-displaymanager)
    - [6.6 Deep Serialization via SessionManager](#66-deep-serialization-via-sessionmanager)
    - [6.7 UndoRedoManager & Parameter History](#67-undoredomanager--parameter-history)
    - [6.8 Fallback Parameter System](#68-fallback-parameter-system)
    - [6.9 Global Tooltip System](#69-global-tooltip-system)
7. [Input Data Format](#7-input-data-format)
8. [GUI Layout Reference](#8-gui-layout-reference)
9. [Tab-by-Tab Technical Guide](#9-tab-by-tab-technical-guide)
    - [Tab 1: Filter](#tab-1-filter)
    - [Tab 2: Detect](#tab-2-detect)
    - [Tab 3: Localize](#tab-3-localize)
    - [Tab 4: Track](#tab-4-track)
    - [Tab 5: Post-Process](#tab-5-post-process)
    - [Tab 6: Render](#tab-6-render)
10. [Specialized Panels & UX Features](#10-specialized-panels--ux-features)
    - [10.1 Visual Adjustments Sidebar](#101-visual-adjustments-sidebar)
    - [10.2 Kalman Trust Balance Panel](#102-kalman-trust-balance-panel)
    - [10.3 Vessel Masking & Interactive Histogram](#103-vessel-masking--interactive-histogram)
    - [10.4 Advanced Parameter Modals](#104-advanced-parameter-modals)
11. [Consolidated Parameter Reference](#11-consolidated-parameter-reference)
12. [Session Management](#12-session-management)
13. [Undo / Redo System](#13-undo--redo-system)
14. [Keyboard Shortcuts](#14-keyboard-shortcuts)
15. [Developer's Guide: Extending the Registry](#15-developers-guide-extending-the-registry)
16. [Performance Optimization](#16-performance-optimization)
17. [Troubleshooting](#17-troubleshooting)
18. [License](#18-license)
19. [Acknowledgments](#19-acknowledgments)

---

## 1. Executive Summary

The **ULM Master GUI v3.0** is a MATLAB-based interactive platform for end-to-end Ultrasound Localization Microscopy (ULM) processing. It was developed to address a fundamental challenge in ULM research: the high sensitivity required to detect and track sub-micron nanobubbles (NB), which have a substantially weaker acoustic backscatter signal compared to standard microbubbles (MB).

Unlike static batch scripts where all parameters are set once and applied globally, the ULM Master GUI provides a **real-time feedback loop**: every parameter adjustment triggers an immediate re-calculation on the cached dataset, visualizing the result instantaneously. This interactive approach is essential for the systematic optimization of the full processing chain — from clutter filter thresholds to sub-pixel localization constraints to tracking gate sizes — and is particularly critical when working with agents operating near the detection limit.

**Core design principles include:** independently tunable processing stages with live visual feedback; a cached SVD decomposition that renders slider adjustments instantaneous; direct export of optimized parameters to the batch processing pipeline (`run_ULM_Analysis_Kidney.m`); and full session save/load for reproducibility across experiments.

**Major features introduced in v3.0** are: a Central Algorithm Registry (`getAlgorithmRegistry`) as the single source of truth for all dropdowns and dispatch; Advanced Parameter Modals exposing every parameter from `setDefaultParams.m` in grouped, annotated modal dialogs; a Global Tooltip System; a Kalman Trust Balance Panel for real-time visualization of model-vs-measurement confidence; a Localization Density Map preview; a persistent Visual Adjustments sidebar; a Robust Fallback System to prevent crashes if `setDefaultParams.m` is absent; and a QC Summary Dialog after localization.

---

## 2. Scientific Background

ULM is a super-resolution imaging technique that reconstructs microvascular architecture by localizing individual gas-filled contrast agents with sub-pixel precision and linking their positions into trajectories over thousands of frames. Resolution is not limited by the acoustic diffraction limit (typically ~λ/2, hundreds of μm for clinical transducers) but by localization precision, which follows:

$$\sigma_{loc} \propto \frac{FWHM}{\sqrt{SNR}}$$

The pipeline consists of five sequential stages:

| Stage | Scientific Challenge | GUI Module |
|-------|---------------------|------------|
| **Clutter Filtering** | Separating bubble signal from tissue background via spatiotemporal decomposition | Tab 1: Filter |
| **Detection** | Identifying candidate bubble peaks above noise in the filtered volume | Tab 2: Detect |
| **Localization** | Estimating sub-pixel bubble centroid using point-spread-function models | Tab 3: Localize |
| **Tracking** | Linking localizations across frames into coherent trajectories | Tab 4: Track |
| **Rendering** | Projecting trajectories onto a high-resolution super-resolved grid | Tabs 5–6 |

This GUI was validated on in vitro gelatin phantoms (100–500 μm channels), in vivo rat brain data (benchmark dataset, Chavignon et al., 2023), and in vivo rat kidney data. Full experimental details are available in the associated thesis (see Citation).

---

## 3. System Requirements

### MATLAB Version

Minimum: MATLAB R2020b. Recommended: MATLAB R2022a or later.

### Required Toolboxes

| Toolbox | Used For |
|---------|----------|
| Image Processing Toolbox | Masking, ROI tools, morphological operations, `adapthisteq`, `imtophat`, `imsharpen` |
| Signal Processing Toolbox | Butterworth filter design, Savitzky-Golay smoothing (`sgolayfilt`) |
| Statistics and Machine Learning Toolbox | K-means clustering in DCC-SVD |
| Optimization Toolbox | Gaussian fitting via `lsqcurvefit` |

### Optional Toolboxes

| Toolbox | Benefit |
|---------|---------|
| Parallel Computing Toolbox | Significant speedup in Gaussian fitting and post-processing via `parfor` |

### Hardware

RAM: 16 GB minimum; 32 GB+ recommended for in vivo datasets. Display: 1920×1080 minimum resolution (GUI is designed for 1600×1000 px). The GUI runs on Windows 10/11 (primary development platform), macOS, and Linux.

---

## 4. Installation & Launch

Download or clone all files into a single directory. Open MATLAB and navigate to the cloned directory. The GUI automatically adds the codebase to the MATLAB path on launch via `addpath(genpath(...))`. To verify toolbox availability, run `ver` in the command window and confirm the toolboxes listed in Section 3 are present.

To launch the GUI:

```matlab
ULM_Master_GUI_v3
```

On launch, the GUI will: attempt to load default parameters via `setDefaultParams.m`; run `ensureAllParamFields()` to fill any missing parameter fields; fall back to safe system defaults via `createFallbackParams()` if `setDefaultParams.m` is missing or fails; display a warning dialog if the fallback was used; open a 1600×1000 px figure window; and apply hover tooltips to all controls via the global tooltip system.

---

## 5. File Structure

```
ulm-codebase/
│
├── [GUI/]
│   ├── ULM_Master_GUI_v3.m       ← Main GUI application (entry point, ~4700 lines)
│   ├── ULM_Constants.m           ← All default values, limits, and validation (centralized)
│   ├── SessionManager.m          ← Save/load session logic with deep serialization
│   ├── UndoRedoManager.m         ← Parameter history management (recursive deep copy)
│   ├── DisplayManager.m          ← Visualization rendering with debounced updates
│   ├── DataHash.m                ← MD5-based change detection for SVD caching
│   └── setDefaultParams.m        ← Experiment-specific parameter configuration
│
├── [core/]
│   ├── ULM_Processor.m           ← Main ULM class
│   ├── [detection/]
│   │   ├── detectBubbles.m       ← Intensity-based regional maxima detection
│   │   ├── detectBubbles_NCC.m   ← Normalized cross-correlation detection
│   │   └── detectBubbles_NP.m    ← Neyman-Pearson hypothesis test detection
│   ├── [localization/]
│   │   ├── localizeRadialSymmetry.m  ← Gradient-based radial symmetry
│   │   ├── fit2DGaussian.m       ← Full 2D Gaussian NLLS fit
│   │   └── fit2DGaussian_Fast.m  ← Vectorized fast Gaussian fitting
│   ├── [tracking/]
│   │   ├── trackNearestNeighbor.m    ← Greedy nearest-neighbor linker
│   │   ├── trackHungarian.m          ← Global Hungarian assignment
│   │   ├── trackKalman.m             ← Standard Kalman filter tracker
│   │   ├── trackKalman_Advanced.m    ← Hierarchical Kalman tracker (HKT)
│   │   ├── calculateCostMatrix.m     ← Smart Cost Matrix (SCM) engine
│   │   └── munkres.m                 ← Munkres/Hungarian solver
│   ├── [filtering/]
│   │   ├── SVD_filter.m          ← Standard SVD clutter filter
│   │   ├── SVD_SSM.m             ← Spatial Similarity Matrix thresholding
│   │   ├── SVD_blockwise.m       ← Block-wise adaptive SVD
│   │   ├── DCC_SVD.m             ← Density Canopy Clustering SVD
│   │   ├── run_SVD_Decomposition.m
│   │   ├── reconstruct_SVD_Signal.m
│   │   └── Butterworth_bandpass_filter.m
│   └── [rendering/]
│       ├── renderHistogram.m     ← Histogram accumulation
│       └── renderGaussian.m      ← Gaussian-blurred density map
│
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
    └── ...
```

---

## 6. Architectural Deep-Dive

### 6.1 Design Patterns & Separation of Concerns

The GUI follows a strict **separation of concerns** with four dedicated manager classes, a centralized algorithm registry, and one monolithic main file that houses all layout builders, callbacks, and the registry itself:

```
ULM_Master_GUI_v3.m          (Layout builder + all callbacks + algorithm registry)
       │
       ├── getAlgorithmRegistry()  (Central registry: filters, trackers, detectors,
       │                            localizers, renderers, smoothers)
       ├── ULM_Constants.m        (Read-only: all defaults, limits, validation logic)
       ├── SessionManager.m       (Save/load: serializes the full app struct)
       ├── UndoRedoManager.m      (History: up to 20 parameter states, recursive deep copy)
       └── DisplayManager.m       (Visualization: frame rendering, overlays, debouncing)
```

The application state is stored in a single `app` struct that is persisted via `guidata(fig, app)` on the main `uifigure` handle. This struct has three primary sub-structures: `app.data` (all raw data, intermediate results, cached SVD matrices, and the `params` struct), `app.ui` (handles to every widget, panel, and axes), and `app.state` (the pipeline state machine, current frame index, playback state, and timer references). All four manager objects are stored directly on `app` as `app.sessionManager`, `app.undoManager`, `app.displayManager`, and `app.constants`.

### 6.2 The Central Algorithm Registry

The `getAlgorithmRegistry()` function is the architectural keystone of v3.0. It returns a struct `reg` with six categories — `reg.filter`, `reg.track`, `reg.detect`, `reg.loc`, `reg.render`, and `reg.smoothing` — each of which is a struct array defining every available algorithm in that category. Each entry carries the following fields:

| Field | Purpose |
|-------|---------|
| `id` | Machine-readable identifier (e.g., `'svd_filter'`, `'Kalman'`). Used as the dropdown value and dispatch key. |
| `display` | Human-readable name shown in dropdown menus (e.g., `'Global SVD'`, `'Hierarchical Kalman'`). |
| `func` | Function name string for dispatch (e.g., `'run_SVD_Decomposition'`, `'trackKalman_Advanced'`). |
| `panel_field` | Name of the UI panel to show/hide when this method is selected (e.g., `'p_svd'`, `'p_dcc'`). |
| `tooltip` | Plain-English description attached to the dropdown via the tooltip system. |
| `isKalman` / `isGaussian` | Boolean flags controlling conditional UI visibility (Kalman settings panel, Gaussian QC panel). |
| `usesHK` / `showsGain` | Flags for the Hierarchical Kalman and Kalman Trust Balance panel. |

All dropdowns are populated from the registry via `registryIds(reg.filter)`, which extracts the `id` fields. Dispatch switches in `runFilter`, `runDetection`, `runLocalization`, and `runTracking` match against these same IDs. UI panel visibility is toggled by `updateFilterOptions`, `updateTrackingOptions`, etc., which iterate the registry to find the entry matching the current dropdown value and then show/hide the associated panel.

> **Tip:** Adding a new algorithm to the GUI requires exactly one step — add a single entry to the appropriate registry array. Dropdowns, dispatch, tooltips, and option panel visibility all derive from the registry automatically. See Section 15 for a step-by-step guide.

### 6.3 Pipeline State Machine & Downstream Invalidation

The GUI state is tracked in `app.state.currentState`, which controls which panels and buttons are enabled. States advance sequentially:

| Stage Code | Meaning | Unlocked Controls |
|:----------:|---------|-------------------|
| -1 | No data loaded | Only Load button active |
| 0 | Data loaded | Filter tab unlocked |
| 1 | Filtered | Detect tab unlocked |
| 2 | Detected | Localize tab unlocked |
| 3 | Localized | Track tab unlocked |
| 4 | Tracked | Post-Process tab unlocked |
| 5 | Post-processed | Render tab unlocked |

The function `manageGUIState(app, stage)` inspects which data fields are populated (e.g., `~isempty(app.data.filteredData)`) and enables or disables the corresponding "Run" buttons accordingly. It uses `matlab.lang.OnOffSwitchState` for safe boolean-to-enable conversion.

**Downstream Invalidation via `clearDownstreamData`:** When a user modifies an early-stage parameter (say, re-running the clutter filter), all subsequent results become stale. The function `clearDownstreamData(data, stage)` enforces data integrity by clearing every field downstream of the given stage. For example, calling `clearDownstreamData(app.data, 1)` after filtering will set `candidateBubbles`, `localizations`, `tracks_raw`, and `tracks_final` to `[]`, forcing the user to re-run detection, localization, tracking, and post-processing. This cascade prevents silent inconsistencies where downstream results do not correspond to the current upstream data.

The invalidation logic uses a simple waterfall:

```
Stage < 1  →  clear filteredData, candidateBubbles, localizations, tracks_raw, tracks_final
Stage < 2  →  clear candidateBubbles, localizations, tracks_raw, tracks_final
Stage < 3  →  clear localizations, tracks_raw, tracks_final
Stage < 4  →  clear tracks_raw, tracks_final
Stage < 5  →  clear tracks_final
```

Every processing callback (e.g., `runFilter`, `runDetection`) calls `clearDownstreamData` after storing its own results, then calls `manageGUIState` to update button states.

### 6.4 SVD Caching via DataHash

SVD decomposition is the most computationally expensive step in the pipeline. For a 122×260×1500 dataset, the Casorati matrix is 31,720×1500, and computing its full SVD can take tens of seconds. However, once computed, the SVD need not be recomputed if only the cutoff range changes — the user can reconstruct any band-pass filtered signal as:

$$X_{filtered} = U \cdot \text{diag}(S_{[c_1:c_2]}) \cdot V^*$$

where $c_1$ and $c_2$ are the cutoff start and end indices. This reconstruction is an $O(N \cdot K)$ matrix multiply (where $K = c_2 - c_1$ is typically small), making it effectively instantaneous.

The caching mechanism uses `DataHash.m`, which computes an MD5 hash of the raw data array using Java's `java.security.MessageDigest`. The hash is stored in `app.data.rawDataHash` as a composite key combining the data hash and the filter method (e.g., `'a3f7c1...2d_svd_filter'`). Before computing a new SVD, the filter callback checks whether `app.data.U` is non-empty and the stored hash matches the current data+method. If so, it skips directly to `reconstruct_SVD_Signal`. For complex (IQ) data, `DataHash` interleaves the real and imaginary parts before hashing to ensure sensitivity to phase changes.

The `DataHash` function also includes a fallback for environments where the JVM is unavailable (e.g., compiled MATLAB): it constructs a fingerprint from the array dimensions, class, first/last element, and sum.

> **Tip:** Slider adjustments to SVD cutoffs after the first filter run are instantaneous because they only call `reconstruct_SVD_Signal` on the cached `U`, `S`, `V` matrices. Only changing the raw data or the filter method triggers a full recomputation.

### 6.5 Debounced Display & the DisplayManager

`DisplayManager` is a `handle` class that encapsulates all visualization logic. Its primary entry point, `displayFrame(app)`, routes rendering to specialized methods based on `app.state.currentState`: `displayRawData` (state 0), `displayFilteredData` (state 1), `displayDetectionResults` (state 2 — overlays red `×` markers), `displayLocalizations` (state 3 — overlays red dots), `displayRawTracks` (state 4 — colored polylines on a mean projection), and `displayProcessedTracks` (state 5 — smoothed polylines filtered by minimum length).

Every displayed image passes through `processImageForDisplay(app, rawImg)`, which applies the following chain from the Visual Adjustments sidebar: absolute value → optional log compression ($20 \cdot \log_{10}(|x| + \epsilon)$) → optional mat2gray normalization → gamma stretch ($x^\gamma$). After image display, `applyAxesDisplaySettings` sets the colormap from the dropdown and either auto-scales or applies manual CLim.

The **debouncing mechanism** protects the UI thread from being overwhelmed by rapid slider movements. The `debouncedDisplay(fig, delay)` method creates a single-shot `timer` with a configurable delay (default: `ULM_Constants.DEBOUNCE_DELAY = 50 ms`). Each new call cancels the previous timer and starts a fresh one. Only when the user pauses slider movement for 50 ms does the actual display refresh fire. This prevents dozens of redundant `imagesc` calls during a single drag gesture.

The `DisplayManager` also pre-computes and caches a `lines(256)` colormap for track rendering to avoid recomputation on every frame.

### 6.6 Deep Serialization via SessionManager

`SessionManager` is a `handle` class responsible for creating complete workspace snapshots and restoring them. Its `createSessionData(app)` method serializes the following into a single struct: version string (`'3.0'`), timestamp, the complete `params` struct, all data arrays (raw data, filtered data, candidate bubbles, localizations, raw tracks, final tracks, mask, SVD components, DCC cluster indices, block-wise diagnostics), and the current GUI state (pipeline stage, frame index, selected tab, slider limits).

The v3.0 `SessionManager` includes a dual-interface design: the original `createSessionData`/`restoreSessionData` pair operates directly on the `app` struct, while the newer `createSession`/`restoreSession` pair provides a flatter API that separates `data` and `params` outputs. The `restoreSession` method reconstructs all derived quantities (color limits via percentile computation, mean projections, etc.) that are not explicitly stored in the session.

Session files use MATLAB's `-v7.3` format (HDF5 backend) to support variables larger than 2 GB. After loading, `ensureAllParamFields()` is called to fill any fields that may be missing in sessions saved with older versions of the codebase, ensuring forward compatibility.

### 6.7 UndoRedoManager & Parameter History

`UndoRedoManager` is a `handle` class that maintains two cell-array stacks — `undoStack` and `redoStack` — each entry containing a timestamped deep copy of the full `params` struct and an operation label. The maximum stack depth is configurable via `ULM_Constants.MAX_UNDO_STATES` (default: 20).

The **deep copy** is critical because MATLAB structs use copy-on-write semantics that can lead to unexpected aliasing when handle objects are involved. The `deepCopy` method recursively traverses the struct tree, copying every nested struct and cell array element by value.

The v3.0 `undo(currentParams)` method pushes the caller-supplied `currentParams` onto the redo stack before popping the previous state from the undo stack, enabling clean round-trip undo/redo. Any new parameter change (via `push`) clears the redo stack, following the standard "branching" undo model.

After undo/redo, `populateGUIFromParams(fig)` is called to synchronize all UI controls to the restored parameter values, using the `setSafe` helper which gracefully handles type mismatches, missing controls, and out-of-range slider values.

### 6.8 Fallback Parameter System

If `setDefaultParams.m` is absent from the MATLAB path or throws an error, the GUI constructs a complete fallback parameter struct via `createFallbackParams()`. This function defines every field the GUI expects across all pipeline stages — IO, processing, acquisition, filter, localization/detection, tracking (including full Kalman and QC sub-structs), rendering, and analysis — with safe, conservative default values (e.g., FPS: 200, pixel size: 0.05 mm, SVD cutoff: [5, 100]). After construction, `ensureAllParamFields()` fills any remaining gaps. A warning dialog informs the user that fallback defaults are in use.

### 6.9 Global Tooltip System

The function `getTooltipDictionary()` returns a struct mapping UI handle names (e.g., `'SVDCutoffStart'`, `'KalmanNoise'`) to plain-English descriptions that include the pipeline stage each parameter belongs to and its effect. At startup, `applyTooltips(app)` iterates every entry in this dictionary and, if a matching handle exists in `app.ui` with a `Tooltip` property, attaches the description. This provides a built-in mini user-guide without leaving the interface.

---

## 7. Input Data Format

### Required Format

The GUI accepts `.mat` files containing any 3D numeric array of shape `[Nz × Nx × Nt]`, where `Nz` is the number of axial pixels, `Nx` is the number of lateral pixels, and `Nt` is the number of frames. The variable name is irrelevant — the loader (`loadRawData`) automatically scans the file for the first 3D numeric field. Both real and complex (IQ) data types are supported; complex data is displayed using `abs()`.

### Physical Calibration

Physical pixel sizes and frame rate are set in the menu bar's **Fundamental Parameters** section:

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `FPS (Hz)` | Acquisition frame rate | 200–1000 Hz |
| `Px X (mm)` | Lateral pixel size | ~0.05–0.1 mm |
| `Px Z (mm)` | Axial pixel size | ~0.05–0.1 mm |

These values propagate to all downstream physical calculations: velocities in mm/s, linking distances in mm, scale bars, and the super-resolution grid.

---

## 8. GUI Layout Reference

The interface is divided into three functional zones:

```
┌─────────────────────────────────────────────────────────────────┐
│  Menu Bar: [Load] [Save] [Undo] [Redo]  [Memory]  FPS: Px:    │
├─────────┬──────────────────────────────┬────────────────────────┤
│ Visual  │                              │  Workflow Tabs:        │
│ Adjust. │                              │  1.Filter 2.Detect ... │
│ Panel   │                              │  ──────────────────    │
│         │   Visualization Canvas       │  Parameter             │
│ Norm.   │   (Central, real-time)       │  Control Panel         │
│ Log     │                              │  (Active tab controls) │
│ Gamma   │                              │                        │
│ Cmap    │                              │  [Advanced... Btn]     │
│ CLim    │                              │  [Run Button]          │
│         │                              │  [Reset to Defaults]   │
├─────────┴──────────────────────────────┴────────────────────────┤
│  Frame Slider ────────────────────────────────── [Frame: ___]   │
│  [Load Data (IQ / imageData)]    [▶ Play / Pause]   Status ●   │
└─────────────────────────────────────────────────────────────────┘
```

### Menu Bar Elements

| Element | Function |
|---------|----------|
| Load Work Session | Restores a previously saved complete GUI state from a `.mat` file |
| Save Work Session | Saves all data, parameters, and results to a `.mat` file (HDF5, `-v7.3`) |
| Undo | Reverts last parameter change (up to 20 levels) |
| Redo | Re-applies an undone change |
| Memory: X GB | Real-time RAM usage monitor (auto-updates every `ULM_Constants.MEMORY_UPDATE_INTERVAL` = 2 s) |
| FPS | Acquisition frame rate — editable in-place, propagates to `params.acq.framerate` |
| Px X / Px Z (mm) | Pixel sizes — editable in-place, propagate to `params.track.pixel_X_size` / `pixel_Z_size` |

### Status Bar

The status lamp and label in the bottom-right indicate the current state: green (Ready / operation complete), red (Processing in progress), blue (Informational, e.g., "Undo complete"), orange (Warning).

---

## 9. Tab-by-Tab Technical Guide

### Tab 1: Filter

**Functional Purpose:** Remove stationary tissue clutter (high-energy, temporally coherent) from the IQ data, retaining the dynamic microbubble signal. Tissue clutter dominates the first singular values of the spatiotemporal Casorati matrix, while noise occupies the last. The bubble signal resides in the intermediate band.

#### Step A: Spatial Crop (Optional, Pre-SVD)

Spatial cropping reduces the Casorati matrix dimensions before SVD, yielding significant speedups. For a 122×260 frame, the Casorati matrix has 31,720 rows. Cropping to a 60×130 ROI reduces this to 7,800 rows — a 4× reduction that translates to approximately 16× SVD speedup because SVD complexity is super-linear in the matrix dimensions.

| Control | Description | Source |
|---------|-------------|--------|
| **Crop Box [x y w h]** | Manual entry of crop rectangle coordinates (pixels) | `app.ui.CropBoxField` |
| **Interactive Crop** | Opens a figure with `drawrectangle` for interactive ROI selection | `runInteractiveCrop` |
| **Load Crop** | Load a previously saved `cropBox.mat` | `loadCrop` |
| **Save Crop** | Export the current crop rectangle for reuse | `saveCrop` |
| **Apply Crop to Data** | Permanently crops `rawData` — destructive within session | `applyCrop` |

> **Warning:** Applying a crop is a destructive operation within the current session. The raw data is replaced by the cropped sub-region, the SVD cache is invalidated, the mask is cleared, and the undo history is reset. A confirmation dialog is shown before proceeding.

#### Step B: Clutter Filter Method

Select from four methods, each populated from `getAlgorithmRegistry().filter`:

**`svd_filter` — Global SVD Filter.** The baseline method. Decomposes the Casorati matrix $X = U \Sigma V^*$ and retains only singular values within the specified range $[c_1, c_2]$:

$$X_{filtered} = \sum_{i=c_1}^{c_2} \sigma_i \cdot u_i \cdot v_i^*$$

| Parameter | Default (`ULM_Constants`) | Code Field | Description |
|-----------|--------------------------|------------|-------------|
| Cutoff Start | 1 (`DEFAULT_SVD_CUTOFF(1)`) | `params.filter.svd_cutoff(1)` | First singular value index to keep. Indices below are discarded as tissue clutter. |
| Cutoff End | 100 (`DEFAULT_SVD_CUTOFF(2)`) | `params.filter.svd_cutoff(2)` | Last singular value index to keep. Indices above are discarded as noise. |

The full SVD is computed once via `run_SVD_Decomposition` and cached. Subsequent slider adjustments call only `reconstruct_SVD_Signal`, which is instantaneous.

**`svd_ssm` — SVD with Spatial Similarity Matrix (SSM).** Implements the method of Baranger et al. (2023). Calculates the Pearson correlation between absolute values of spatial singular vectors ($U$) to identify the transition boundary between tissue-dominated and bubble-dominated components automatically. Uses the same Cutoff Start/End sliders as the standard SVD but provides a statistical guide for setting the upper cutoff.

**`dcc_svd` — Density Canopy Clustering SVD.** Implements the DCC-SVD method of Han et al. (2024). Each singular component is characterized by a 3D feature vector: normalized log-energy, power-weighted central temporal frequency, and spatial correlation to the mean spatial vector. K-means clustering (seeded by density canopy centers) partitions components into "Tissue", "Blood", and "Noise" clusters. Interactive sliders on the DCC panel allow manual fine-tuning of cluster boundaries after the initial automatic classification. The GUI caches the DCC-computed SVD components (`U`, `S`, `V`) and the cluster indices (`tissue_indices`, `blood_indices`, `noise_indices`) so that boundary adjustments reconstruct the signal instantly.

| Parameter | Description |
|-----------|-------------|
| Tissue Start/End (%) | Percentage range of the tissue cluster |
| Blood Start/End (%) | Percentage range of the blood cluster |
| Noise Start/End (%) | Percentage range of the noise cluster |

**`svd_blockwise` — Block-Wise Adaptive SVD.** Implements the method of Song et al. (2017). Divides the image into overlapping spatial blocks and applies independent adaptive SVD thresholding to each block, accounting for spatially varying clutter characteristics. The dedicated parameter panel appears when this method is selected:

| Parameter | Default | Code Field | Description |
|-----------|---------|------------|-------------|
| Threshold Method | `DopplerGradient` | `params.filter.blockwise.threshold_method` | Strategy for per-block cutoff: `DopplerGradient`, `SSM`, `Hybrid`, or `Manual` |
| Manual Cutoff [Lo Hi] | [10 200] | `params.filter.blockwise.manual_cutoff` | Only visible when Method = `Manual` |
| Block Size (mm) | 4.0 | `params.filter.blockwise.block_size_mm` | Spatial extent of each processing block |
| Overlap (%) | 75 | `params.filter.blockwise.overlap_pct` | Block overlap. 75% for exploratory, 93.75% for publication-quality |
| MP Deviation (σ) | 2.0 | `params.filter.blockwise.mp_deviation_sigma` | Marchenko-Pastur sensitivity for the high cutoff |
| Gradient Inflection (%) | 0.10 | `params.filter.blockwise.gradient_pct` | Sensitivity of inflection detection for Cutoff 1A |
| Tissue Freq Thr (Hz) | -1 (auto) | `params.filter.blockwise.tissue_freq_hz` | Set to -1 for automatic: `max(5, min(20, framerate/50))` |
| Min Blood Comps | 3 | `params.filter.blockwise.min_blood_comps` | Floor constraint on blood components per block |
| Max Tissue Fraction | 0.60 | `params.filter.blockwise.max_tissue_frac` | Ceiling on the fraction classified as tissue |
| Plot threshold maps | Off | `params.filter.blockwise.plot_maps` | Display spatial threshold maps after filtering |

#### Step C: Butterworth Bandpass Filter (Optional)

An auxiliary temporal frequency filter applied after SVD to further isolate the bubble signal within a specific frequency band.

| Parameter | Default | Code Field | Description |
|-----------|---------|------------|-------------|
| Enable | Off | `params.filter.enable_butterworth` | Master switch |
| Cutoff [Low, High] (Hz) | [10, 100] | `params.filter.butter_cutoff` | Passband frequencies |
| Order | 4 | `params.filter.butter_order` | Filter order; higher = sharper rolloff |

Validation is performed via `ULM_Constants.isValidButterCutoff(cutoff, framerate)`, which ensures the upper cutoff is below the Nyquist frequency.

#### Step D: Spatial Filter (Optional)

A per-frame spatial convolution filter for additional noise reduction. Controls are shown and hidden dynamically depending on the selected method:

| Parameter | Options | Default | Description |
|-----------|---------|---------|-------------|
| Method | `None`, `Gaussian`, `Median`, `DoG`, `Top-Hat` | `Gaussian` | Type of spatial filter |
| Kernel Size (px) | — | 3 | Filter kernel size (must be odd). Visible for Gaussian, Median, Top-Hat. |
| Sigma 1 | — | 1.0 | Primary Gaussian sigma. Visible for Gaussian and DoG. |
| Sigma 2 | — | 2.0 | Secondary sigma. Visible for DoG only. |

The `updateSpatialOptions` callback dynamically shows/hides the kernel, sigma1, and sigma2 controls based on the selected method.

---

### Tab 2: Detect

**Functional Purpose:** Identify candidate bubble locations in the filtered frames as integer-pixel regional maxima, constrained by an optional vessel mask to reduce false positives and computation.

#### Step A: ROI Masking (Vessel Map)

The upper half of Tab 2 provides both manual mask loading and algorithmic vessel masking.

**Manual mask controls:** Load Mask (loads a binary `.mat` file), Create New Mask (commits the current ROI preview as the active detection mask), Reset Mask (clears the loaded mask), Save Mask (exports to `.mat`), and a status indicator.

**Algorithmic Auto-Mask (Flow Probability Mask):** Generates a mask from the temporal average of the SVD-filtered sequence. The temporal mean image is computed with square-root compression ($\sqrt{\overline{|x|}}$), then processed through an enhancement pipeline with real-time histogram feedback:

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| Enhancement Method | `None` | `None`, `CLAHE`, `Top-Hat`, `Sharpen` | Image enhancement algorithm |
| Enhancement Strength | 0.5 | 0–1 | Enhancement intensity |
| Gamma | 1.0 | 0.1–3.0 | Contrast adjustment ($I' = I^\gamma$) |
| Threshold | 0.0 | 0–1.1 | Binary threshold for mask generation |

The enhancement methods operate as follows. **CLAHE** (Contrast Limited Adaptive Histogram Equalization) uses `adapthisteq` with the clip limit controlled by the Strength slider: `clipLim = 0.001 + strength × 0.04`, and 8×8 tiles. **Top-Hat** filtering uses `imtophat` with a disk structural element whose radius is `1 + round(strength × 10)`, extracting vessel-like structures. **Sharpen** uses `imsharpen` with `Amount = strength × 2`.

The histogram panel (`axHist`) displays the intensity distribution of the enhanced image in log scale, with a red vertical line at the current threshold value. As the user adjusts any parameter, the enhancement is re-applied in real-time via `applyVesselEnhancement`, the histogram updates, and the mask overlay on the main canvas refreshes.

#### Step B: Detection Parameters

| Parameter | Default (`ULM_Constants`) | Range | Code Field | Description |
|-----------|--------------------------|-------|------------|-------------|
| Detection Method | `Intensity` | `Intensity`, `NP`, `NCC` | `params.loc.DetectMethod` | Algorithm for finding candidate peaks |
| Intensity Threshold | 0.3 (`DEFAULT_DETECTION_THRESHOLD`) | 0.01–1.0 | `params.loc.detection_threshold` | Normalized intensity threshold |
| Max Bubbles/Frame | 2000 (`DEFAULT_MAX_BUBBLES_PER_FRAME`) | 1–5000 | `params.loc.max_bubbles_per_frame` | Maximum candidates per frame |
| PSF FWHM [x z] (px) | [1.5, 1.5] (`DEFAULT_FWHM`) | — | `params.loc.fwhm` | Estimated PSF size in pixels |
| Peak Contrast (h) | 0 | 0–0.5 | `params.loc.h_contrast` | H-maxima contrast for dense-field mode |

**Detection Methods:**

**Intensity** (`detectBubbles.m`): Regional maxima above a normalized intensity threshold. Fast and robust. The threshold is applied as a fraction of the frame's maximum intensity.

**Neyman-Pearson** (`detectBubbles_NP.m`): A hypothesis-test detector with controlled false alarm rate. The null hypothesis is that a pixel belongs to background noise. The key parameter `NP_alpha0` (default: 0.01) sets the false alarm probability $P_{FA} = \alpha_0$, from which the detection threshold is derived from the noise statistics.

**Normalized Cross-Correlation** (`detectBubbles_NCC.m`): Template matching against a reference PSF kernel. If no PSF template (`MB_image`) is loaded, the GUI offers to auto-generate one from the current FWHM setting: $\text{PSF}(x,z) = \exp\left(-\frac{x^2}{2\sigma_x^2} - \frac{z^2}{2\sigma_z^2}\right)$ where $\sigma = \text{FWHM}/2.355$. The key parameter `NCC_tau` (default: 0.7) sets the minimum correlation coefficient for peak acceptance.

> **Tip (Saturation Check):** If the status bar reports a detection count exactly equal to `Max Bubbles/Frame`, the detector is saturated and valid bubbles are being missed. Increase the limit until the detected count drops below it.

---

### Tab 3: Localize

**Functional Purpose:** Refine integer-pixel candidate positions to sub-pixel precision using the PSF geometry.

#### Localization Methods

| Method | Algorithm | Speed | Precision | Registry Flag |
|--------|-----------|-------|-----------|---------------|
| `radial` | Gradient-based radial symmetry (Parthasarathy, 2012) | Fast (non-iterative) | High | `isGaussian: false` |
| `gaussian_fit` | 2D Gaussian NLLS (Levenberg-Marquardt) | Slow (iterative) | Highest | `isGaussian: true` |
| `gaussian_fit_fast` | Vectorized/linearized Gaussian fitting | Medium | High | `isGaussian: true` |

**Radial Symmetry Localization.** The sub-pixel center $(\hat{x}, \hat{y})$ is the weighted least-squares intersection of gradient lines emanating from each pixel in the ROI:

$$\text{weights} = |\nabla I|^2$$

Pixels on the PSF slopes carry the most information (highest gradient magnitude), naturally weighting the center estimate toward the highest-SNR pixels. This is a closed-form, non-iterative solution requiring no initial guess. Key internal thresholds: `min_gradient_squared` ($10^{-6}$, rejects flat regions) and `min_determinant` ($10^{-6}$, rejects numerically unstable fits).

**Gaussian Fit Localization.** Fits the parametric model:

$$I(x,y) = A \cdot \exp\left(-\frac{(x-x_0)^2 + (y-y_0)^2}{2\sigma^2}\right) + C$$

to each ROI using nonlinear least-squares (`lsqcurvefit`). Reports goodness-of-fit ($R^2$) for each localization.

#### Localization Quality Control (QC)

A multi-layer filter removes candidates that do not conform to the expected PSF shape:

| QC Check | Default | Code Field | Description |
|----------|---------|------------|-------------|
| Divergence Check | On | `params.loc.enable_divergence_check` | Rejects solutions shifting more than `max_shift_factor × FWHM/2` from the coarse peak |
| Max Shift Factor | 2.0 | `params.loc.qc_max_shift_factor` | Multiplier on allowable shift |
| ROI Maxima Check | On | `params.loc.enable_roi_maxima_check` | Rejects ROIs with multiple intensity peaks (overlapping bubbles) |
| Min R² (Gaussian only) | 0.3 | `params.loc.min_r_squared` | Minimum goodness-of-fit for Gaussian fitting |
| Box Radius (Gaussian) | 3 px | `params.loc.gauss_box_radius` | Half-size of the fitting ROI |

#### Advanced Localization Parameters (Modal)

Click **"Advanced Detection / PSF Parameters..."** to access:

| Parameter | Default | Description |
|-----------|---------|-------------|
| PSF Type | `Gaussian` | Template type for NCC detection |
| PSF Size [x z] | [5, 5] | Template dimensions |
| Max ROI Maxima | 3 | Maximum local maxima allowed in an ROI before rejection |
| Min |∇|² for fit | 1e-6 | Minimum gradient magnitude for radial symmetry |
| Min Hessian determinant | 1e-6 | Minimum matrix determinant for radial symmetry linear system |

#### Localization Density Map & QC Summary

Click **"Show Localization Density Map"** to render a pre-tracking preview of localization density. The map uses the current upsampling factor and displays accumulated counts per super-resolved pixel with power-law compression ($\gamma = 1/3$) and a `hot` colormap. After localization completes, a QC Summary Dialog displays rejection counts per criterion, total pass/fail rates, and the captured console output in a monospaced text area.

> **Tip (Yield Monitoring):** A drop from 141 detections to 107 localizations (24% rejection) is typical and healthy. A rejection rate exceeding 60% suggests the FWHM or QC parameters need adjustment.

---

### Tab 4: Track

**Functional Purpose:** Link sub-pixel localizations across frames into continuous trajectories using predictive motion models and assignment solvers.

#### Tracking Algorithms

| Algorithm | Registry ID | File | Motion Model | Prediction |
|-----------|------------|------|-------------|------------|
| Kalman (KT) | `Kalman` | `trackKalman.m` | Constant-velocity state $[x, y, v_x, v_y]$ | Yes |
| Hierarchical Kalman (HKT) | `Kalman_Advanced` | `trackKalman_Advanced.m` | Multi-pass velocity bands | Yes |
| Hungarian (HT) | `Hungarian` | `trackHungarian.m` | None (Munkres assignment) | No |
| Nearest Neighbor (NN) | `nn` | `trackNearestNeighbor.m` | None (greedy) | No |

#### Core Tracking Parameters

| Parameter | Default (`ULM_Constants`) | Range | Code Field | Description |
|-----------|--------------------------|-------|------------|-------------|
| Max Linking Distance (px) | 2.0 (`DEFAULT_LINKING_DISTANCE`) | 0.1–10 | `params.track.max_linking_distance` | Maximum Euclidean distance for linking |
| Max Gap Closing (frames) | 2 (`DEFAULT_GAP_CLOSING`) | 0–10 | `params.track.max_gap_closing_frames` | Consecutive missing frames before termination |
| Min Track Length | 3 (`DEFAULT_TRACK_LENGTH`) | 2–20 | `params.track.min_track_length` | Minimum localizations to retain a track |

#### Kalman Filter Parameters

| Parameter | Default | Code Field | Description |
|-----------|---------|------------|-------------|
| Motion Model | `ConstantVelocity` | `params.track.kalman.motion_model` | State vector model |
| Process Noise | 10 | `params.track.kalman.process_noise` | Model flexibility. Low = stiff, high = follows measurements |
| Assignment Method | `hungarian` | `params.track.kalman.assignment_method` | Inner solver: `hungarian` (optimal) or `nn` (fast) |

#### Smart Cost Matrix (SCM)

When **Use Advanced Cost Matrix** is enabled, the linking cost incorporates direction and brightness consistency:

$$C_{total}(i,j) = C_{dist} \cdot (1 + W_{slope} \cdot P_{angle}) \cdot (1 + W_{int} \cdot P_{intensity})$$

where $C_{dist}$ is the base spatial distance (Euclidean or Kalman-predicted), $P_{angle} = \max(0, \theta_{ij} - \theta_{gate})$ is the directional penalty for turns beyond the safety gate, and $P_{intensity} = |I_{current} - \bar{I}_{track}| / (\bar{I}_{track} + \epsilon)$ is the brightness consistency penalty.

The Advanced Kalman modal exposes:

| Parameter | Default | Code Field | Description |
|-----------|---------|------------|-------------|
| Direction penalty weight | 2 | `params.track.kalman.direction_penalty_weight` | Weight on direction deviation (0 = disabled) |
| Angle penalty slope | 0.3 | `params.track.kalman.angle_penalty_slope` | Slope of directional cost above the gate angle |
| Brightness penalty weight | 2 | `params.track.kalman.brightness_penalty_weight` | Weight on amplitude mismatch |
| Max angle change (°) | 70 | `params.track.kalman.max_angle_change_deg` | Hard ceiling on direction change |
| Gating angle change (°) | 90 | `params.track.kalman.gating_max_angle_change_deg` | Soft gate for pre-filtering |
| Direction history points | 4 | `params.track.kalman.direction_history_points` | Past positions used for direction estimation |

#### Hierarchical Kalman Tracker (HKT)

The HKT decomposes tracking into $N$ velocity levels, processing slow bubbles first and subtracting their localizations before tracking faster ones. Parameters:

| Parameter | Default | Code Field | Description |
|-----------|---------|------------|-------------|
| HK alpha | 0.01 | `params.track.kalman.hk_alpha` | Process noise = α × v_max_level |
| HK beta | 0.025 | `params.track.kalman.hk_beta` | Measurement noise = β / 2^(level−1) |
| Max Velocity (mm/s) | 20 | `params.track.kalman.hk_v_max` | Global upper velocity limit |
| Num Levels | 5 | `params.track.kalman.hk_num_levels` | Number of velocity bands |
| Spacing Power | 1.0 | `params.track.kalman.hk_spacing_power` | 1.0 = linear; >1.0 = more levels at low velocities |
| Enable Overlap | On | `params.track.kalman.hk_enable_overlap` | Overlap band between velocity levels |
| Overlap Width (mm/s) | 2.0 | `params.track.kalman.hk_overlap_mm_s` | Width of the overlap band |
| Forward-Backward | On | `params.track.kalman.hk_forward_backward` | Dual-pass: forward then backward |

#### Post-Tracking Quality Control

| QC Constraint | Default | Code Field | Description |
|--------------|---------|------------|-------------|
| Direction Constraint | On | `params.track.qc.enable_direction_constraint` | Rejects tracks with excessive turns |
| Max Angle (°) | 90 (`DEFAULT_MAX_ANGLE_CHANGE`) | `params.track.qc.max_angle_change_deg` | Maximum allowed direction change |
| Acceleration Constraint | On | `params.track.qc.enable_acceleration_constraint` | Rejects unrealistic velocity jumps |
| Factor C | 3.0 (`DEFAULT_ACCELERATION_C_FACTOR`) | `params.track.qc.acceleration_C_factor` | Acceleration tolerance factor |
| VD (Jitter) Constraint | Off | `params.track.qc.enable_vd_constraint` | Rejects tracks with high tortuosity |
| Max VD Ratio | 2.0 (`DEFAULT_VD_RATIO`) | `params.track.qc.max_vd_ratio` | Maximum path-length / displacement ratio |

---

### Tab 5: Post-Process

**Functional Purpose:** Smooth discrete trajectory points and interpolate onto a sub-pixel grid before rendering.

| Parameter | Default | Range | Code Field | Description |
|-----------|---------|-------|------------|-------------|
| Enable Smoothing | On | — | `params.track.enable_postprocessing` | Master switch |
| Window Size | 5 | 3–21 (odd) | `params.track.smoothing_factor` | Smoothing window width |
| Min Length (display) | — | 2–20 | `params.track.display_min_length` | Minimum track length for display (live, no re-run needed) |

The default smoothing method is Savitzky-Golay (polynomial order 3, configured via `ULM_Constants.SGOLAY_POLY_ORDER`). If the track is shorter than the minimum window for SG fitting (poly_order + 2), a `movmean` fallback is used. After smoothing, the track is interpolated at sub-pixel resolution using `interp1` with the method configured in `params.render.interpolation_method` (default: `pchip`), with `fillmissing` to handle edge cases. Velocities are computed as:

$$v_i = \frac{\sqrt{(\Delta x_i \cdot p_x)^2 + (\Delta z_i \cdot p_z)^2}}{\Delta t_i}$$

where $p_x, p_z$ are the physical pixel sizes and $\Delta t_i$ is the inter-frame time interval.

If the Parallel Computing Toolbox is available, post-processing uses `parfor` over tracks.

#### Smoothing Methods (via Advanced Render modal)

| Method | Shape Preservation | Noise Suppression | Speed |
|--------|-------------------|------------------|-------|
| `sgolay` (Savitzky-Golay) | Excellent | Very Good | Fast |
| `gaussian` | Very Good | Excellent | Fast |
| `movmean` (Moving Average) | Good | Good | Fastest |
| `rloess` (Robust Loess) | Very Good | Good | Slow |

---

### Tab 6: Render

**Functional Purpose:** Project all trajectories onto a high-resolution grid to generate super-resolution images.

#### Rendering Parameters

| Parameter | Default (`ULM_Constants`) | Range | Code Field | Description |
|-----------|--------------------------|-------|------------|-------------|
| Upsampling Factor | 3 (`DEFAULT_UPSAMPLING_FACTOR`) | 1–10 | `params.render.upsampling_factor` | Resolution multiplier relative to native pixel |
| Render Method | `histogram` (`DEFAULT_RENDER_METHOD`) | `histogram`, `gaussian` | `params.render.method` | Accumulation strategy |

#### Advanced Rendering Parameters (Modal)

| Parameter | Default | Code Field | Description |
|-----------|---------|------------|-------------|
| Smoothing method | `sgolay` | `params.track.smoothing_method` | Algorithm for smoothing position traces |
| Interpolation method | `spline` | `params.render.interpolation_method` | Sub-step interpolation: `spline`, `pchip`, `linear`, `makima` |
| Gaussian sigma (px) | 0.3 | `params.render.gaussian_sigma` | Spread of each localization when using Gaussian splatting |
| Interpolation step | 0.5 (`DEFAULT_INTERPOLATION_STEP`) | `params.render.interpolation_step` | Sub-pixel spacing between interpolated points |
| Tortuosity bin step | 0.05 | `params.analysis.tortuosity_bins` | Bin width for tortuosity histogram |
| Velocity histogram bins | 60 | `params.analysis.velocity_hist_num_bins` | Number of bins in velocity histogram |
| Density grid (mm) | 0.5 | `params.analysis.density_grid_size_mm` | Cell size of density map for statistics |

#### Generated Outputs

Clicking **"Generate & Display Final Images (New Windows)"** produces four super-resolution maps simultaneously:

| Map | Algorithm | Description |
|-----|-----------|-------------|
| **Density Map** | `accumarray` with cube-root compression ($\gamma = 1/3$) | Accumulates bubble counts per super-resolved pixel. Power-law compression simultaneously visualizes high-flux arteries and sparse capillaries. |
| **Raw Velocity Map** | Per-pixel mean of instantaneous velocities | Arithmetic mean of $|p_t - p_{t-1}|/\Delta t$ at each pixel. Unbiased statistical representation. |
| **Filtered Velocity Map** | `imgaussfilt` with σ = 0.6 | Gaussian-smoothed version of the raw velocity map. Bridges discrete sampling gaps; suppresses outliers. |
| **Combined Fusion Map** | HSV encoding | Hue = velocity (jet colormap, blue→red), Value = density (power-law compressed). Correlates anatomy with hemodynamics in a single image. |

All maps include automatic 1 mm scale bar annotation in the lower-right corner, with font size and line width scaled to the upsampling factor.

---

## 10. Specialized Panels & UX Features

### 10.1 Visual Adjustments Sidebar

The persistent left-side panel provides display controls that apply to every processing stage. These controls affect only the visualization pipeline within `DisplayManager.processImageForDisplay`, not the underlying data:

| Control | Range | Default | Effect on Rendering Pipeline |
|---------|-------|---------|------------------------------|
| Normalize (mat2gray) | On/Off | On | Normalizes frame to [0, 1] via `mat2gray` |
| Log Compression | On/Off | Off | Applies $20 \cdot \log_{10}(|x| + \epsilon)$ for high dynamic range data |
| Gamma (Stretch) | 0.1–5.0 | 1.0 | Applies $x^\gamma$ — values < 1 brighten dark regions, > 1 darken |
| Colormap | `gray`, `hot`, `jet`, `parula` | `gray` | Sets `colormap(ax, ...)` on the visualization axes |
| Auto CLim | On/Off | On | Automatically computes color limits per frame |
| CLim [Min, Max] | — | [0, 1] | Manual override when Auto CLim is off |

### 10.2 Kalman Trust Balance Panel

This panel (visible on Tab 4 when a Kalman tracker is selected) shows the theoretical split between trusting the motion model vs. trusting the raw localizations, computed from the current noise parameters — no tracking run is needed. The panel updates live as the user adjusts Process Noise, FWHM, or HK alpha/beta.

The computation depends on the selected tracker:

**Standard Kalman:**

$$K = \frac{Q}{Q + R}, \quad Q = \text{process\_noise}, \quad R = \left(\frac{\overline{FWHM}}{2.355}\right)^2$$

**Hierarchical Kalman:**

$$K = \frac{\alpha \cdot v_{max}}{\alpha \cdot v_{max} + \beta}$$

Where $K \to 0$ means the tracker trusts the motion model (smooth, predictive), and $K \to 1$ means it trusts the raw localizations (follows data closely).

The panel displays a **colored split bar** (green = Model, blue = Localizations) with dynamically updated percentages and column widths. The formula text below the bar shows the exact computation with the current parameter values. The color saturation of each bar half intensifies as the balance becomes more extreme, providing an immediate visual cue.

**Practical guidance:** For phantom data with known straight channels, a low $K$ (~0.1–0.3) is appropriate — the model can confidently predict linear motion. For in vivo data with tortuous capillaries, a higher $K$ (~0.5–0.7) allows the tracker to follow the actual vessel geometry more closely. If $K > 0.9$, the Kalman filter provides almost no prediction benefit over a simple nearest-neighbor linker.

### 10.3 Vessel Masking & Interactive Histogram

The ROI panel on Tab 2 provides an integrated vessel masking workflow. When the user navigates to the Detect tab with filtered data available, `prepareROITab` computes a temporal mean image with square-root compression as the base vessel map. The enhancement pipeline (CLAHE / Top-Hat / Sharpen) processes this base image through the user's selected method, applies gamma contrast, and normalizes to [0, 1].

The **histogram panel** (`axHist`) displays the intensity distribution of the enhanced vessel map in log scale. A red vertical threshold line (`xline` with tag `'threshLine'`) moves in real-time as the threshold slider changes — the histogram itself is not redrawn on threshold changes; only the `xline` value is updated, providing smooth interaction. Full histogram redraws occur only when the enhancement method, strength, or gamma changes. Pixels above the threshold are marked as vessel mask (`app.data.mask = vesselMap >= threshold`).

### 10.4 Advanced Parameter Modals

Three "Advanced..." buttons open modal dialogs that expose the full parameter set in grouped, annotated panels. Each dialog follows a consistent visual grammar: grouped `uipanel` elements titled with the pipeline stage, bold parameter labels, numeric edit fields, and italic one-line descriptors. Save and Cancel buttons persist or discard changes. All three modals call `clearDownstreamData` on save to invalidate affected results.

**Advanced Detection / Localization Parameters:** Exposes PSF template configuration (type, size, file path), QC thresholds (max ROI maxima, minimum gradient squared, minimum Hessian determinant), and method-specific controls.

**Advanced Kalman Parameters:** Organized into three sections — Cost Matrix Weights (direction penalty, angle penalty slope, brightness penalty), Angle Gating (max angle change, gating angle, direction history points), and HK Noise Scaling (alpha, beta). The HK alpha and beta fields are wired to `updateKalmanGainSummary` via `ValueChangedFcn`, providing live Trust Balance updates while the dialog is open.

**Advanced Rendering / Smoothing / Analysis Parameters:** Exposes smoothing method, interpolation method, Gaussian sigma, interpolation step, tortuosity bin step, velocity histogram bins, and density grid cell size.

---

## 11. Consolidated Parameter Reference

The following table presents the authoritative default values from `ULM_Constants.m` and `createFallbackParams()`:

| Category | Parameter | Constant Name | Default | Range |
|----------|-----------|--------------|---------|-------|
| **Acquisition** | Frame Rate (Hz) | `DEFAULT_FRAMERATE` | 1000 | — |
| **Acquisition** | Pixel Size X (mm) | `DEFAULT_PIXEL_SIZE_X` | 0.1 | — |
| **Acquisition** | Pixel Size Z (mm) | `DEFAULT_PIXEL_SIZE_Z` | 0.1 | — |
| **Filter** | SVD Cutoff | `DEFAULT_SVD_CUTOFF` | [1, 100] | [1, Nt] |
| **Filter** | Butterworth Cutoff (Hz) | `DEFAULT_BUTTER_CUTOFF` | [10, 100] | (0, framerate/2) |
| **Filter** | Butterworth Order | `DEFAULT_BUTTER_ORDER` | 4 | 1–8 |
| **Detection** | Intensity Threshold | `DEFAULT_DETECTION_THRESHOLD` | 0.3 | 0.01–1.0 |
| **Detection** | Max Bubbles/Frame | `DEFAULT_MAX_BUBBLES_PER_FRAME` | 2000 | 1–5000 |
| **Localization** | FWHM [X,Z] (px) | `DEFAULT_FWHM` | [1.5, 1.5] | — |
| **Localization** | Gauss Box Radius (px) | `DEFAULT_GAUSS_BOX_RADIUS` | 3 | 2–10 |
| **Localization** | Max Shift Factor | `DEFAULT_QC_MAX_SHIFT_FACTOR` | 2.0 | — |
| **Tracking** | Max Linking Distance (px) | `DEFAULT_LINKING_DISTANCE` | 2.0 | 0.1–10 |
| **Tracking** | Max Gap Closing (frames) | `DEFAULT_GAP_CLOSING` | 2 | 0–10 |
| **Tracking** | Min Track Length | `DEFAULT_TRACK_LENGTH` | 3 | 2–20 |
| **Tracking** | Kalman Process Noise | `DEFAULT_KALMAN_PROCESS_NOISE` | 10 | — |
| **Tracking** | Kalman Model | `DEFAULT_KALMAN_MODEL` | `ConstantVelocity` | — |
| **Tracking** | Assignment Method | `DEFAULT_ASSIGNMENT_METHOD` | `hungarian` | — |
| **QC** | Max Angle Change (°) | `DEFAULT_MAX_ANGLE_CHANGE` | 90 | — |
| **QC** | Acceleration C Factor | `DEFAULT_ACCELERATION_C_FACTOR` | 3.0 | — |
| **QC** | VD Ratio | `DEFAULT_VD_RATIO` | 2.0 | — |
| **Post-Process** | Smoothing Window | `DEFAULT_SMOOTHING_WINDOW` | 5 | 3–21 (odd) |
| **Post-Process** | Interpolation Step | `DEFAULT_INTERPOLATION_STEP` | 0.5 | — |
| **Render** | Upsampling Factor | `DEFAULT_UPSAMPLING_FACTOR` | 3 | 1–10 |
| **Render** | Method | `DEFAULT_RENDER_METHOD` | `histogram` | — |
| **ROI** | Gamma | `DEFAULT_ROI_GAMMA` | 1.0 | 0.1–3.0 |
| **ROI** | Threshold | `DEFAULT_ROI_THRESHOLD` | 0.0 | — |
| **UI** | Debounce Delay (s) | `DEBOUNCE_DELAY` | 0.05 | — |
| **UI** | Memory Update (s) | `MEMORY_UPDATE_INTERVAL` | 2.0 | — |
| **UI** | Playback FPS | `PLAYBACK_FPS` | 20 | — |
| **Undo** | Max History States | `MAX_UNDO_STATES` | 20 | — |
| **Performance** | Large Data Threshold | `LARGE_DATA_THRESHOLD` | 1e8 | — |
| **Scale Bar** | Length (mm) | `SCALE_BAR_LENGTH_MM` | 1.0 | — |

---

## 12. Session Management

The session system allows complete workspace preservation and reproducibility.

**Saving a Session:** Click **Save Work Session** in the menu bar. A `.mat` file is created (with `-v7.3` / HDF5 backend for large files) containing: loaded raw data, all processed intermediate results (filtered data, detections, localizations, tracks), all current parameter values, mask and crop settings, GUI state (current tab, frame index), and DCC/block-wise diagnostic data.

**Loading a Session:** Click **Load Work Session**. All data and parameters are restored, `ensureAllParamFields()` fills any missing fields from older sessions, `populateGUIFromParams` synchronizes all UI controls, the correct pipeline stage is inferred from which data fields are populated, and the display updates automatically.

**Session Info (without full load):** From the MATLAB command window:

```matlab
sm = SessionManager();
info = sm.getSessionInfo('my_session.mat');
```

Returns metadata (version, date, processing state, data size, track count) without loading the full dataset.

---

## 13. Undo / Redo System

The undo system tracks parameter changes with named checkpoints:

| Action | Control | Notes |
|--------|---------|-------|
| **Undo** | Undo button or **Ctrl+Z** | Reverts last parameter change |
| **Redo** | Redo button or **Ctrl+Y** | Re-applies undone change |
| **History** | `app.undoManager.displayHistory()` | Prints full undo/redo chain to console |

Up to 20 levels of undo are maintained (configurable via `ULM_Constants.MAX_UNDO_STATES`). Undo/Redo only tracks parameter changes, not processing results — after undoing, click the relevant "Run" button to re-process with the restored parameters. The Undo and Redo buttons are automatically greyed out when no history is available. Applying a spatial crop clears the entire undo history. Making any new parameter change clears the redo stack (standard branching model).

---

## 14. Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| **Ctrl+Z** | Undo last parameter change |
| **Ctrl+Y** | Redo |
| **Spacebar** | Play / Pause frame playback |
| **← Arrow** | Previous frame (when frame slider is focused) |
| **→ Arrow** | Next frame (when frame slider is focused) |

---

## 15. Developer's Guide: Extending the Registry

To add a new algorithm (e.g., a new tracking method called `myTracker`), follow these steps:

**Step 1.** Implement the function `trackMyTracker.m` in `core/tracking/` with the standard signature:

```matlab
function tracks = trackMyTracker(localizations, params, indentPrefix)
    % Your tracking logic here
end
```

**Step 2.** Open `ULM_Master_GUI_v3.m` and locate `getAlgorithmRegistry()`. Add one entry to `reg.track`:

```matlab
reg.track = struct( ...
    'id',         {'Kalman', 'Kalman_Advanced', 'Hungarian', 'nn', 'myTracker'}, ...
    'display',    {'Kalman', 'Hierarchical Kalman', 'Hungarian', 'Nearest Neighbor', 'My Tracker'}, ...
    'func',       {'trackKalman','trackKalman_Advanced','trackHungarian','trackNearestNeighbor','trackMyTracker'}, ...
    'isKalman',   {true, true, false, false, false}, ...
    'showsGain',  {true, true, false, false, false}, ...
    'usesHK',     {false, true, false, false, false}, ...
    'tooltip',    { ...
        'Classic Kalman filter...', ...
        'Hierarchical Kalman...', ...
        'Linear-assignment solver...', ...
        'Greedy nearest-neighbor...', ...
        'My custom tracker description.' ...
    });
```

**Step 3.** Add a `case` to the dispatch switch in `runTracking`:

```matlab
case 'mytracker'
    tracks = trackMyTracker(app.data.localizations, params, '  ');
```

**That's it.** The dropdown will now include "My Tracker", the tooltip will appear on hover, and the existing option visibility logic will handle any panel that needs showing/hiding based on the `isKalman`, `usesHK`, and `showsGain` flags.

---

## 16. Performance Optimization

### For Large Datasets (>100 MB)

**Spatial crop first.** This reduces SVD matrix size quadratically. A 2× crop in each dimension yields a 4× reduction in the Casorati matrix rows, giving approximately 16× SVD speedup. Use "Interactive Crop" then "Apply Crop to Data" in Tab 1.

**Use the cached SVD.** After the first filter run, slider adjustments to SVD cutoffs are instantaneous because they only call `reconstruct_SVD_Signal` on the cached `U`, `S`, `V` matrices. Avoid clicking "Run Filter" again unless the raw data or filter method changes.

**Enable Parallel Computing Toolbox.** The post-processing callback uses `parfor` when a parallel pool is detected:

```matlab
parpool('local', 4);   % Start 4 parallel workers
```

**Use `histogram` render mode.** The `gaussian` render mode is significantly slower for large track sets because it splatts a Gaussian kernel for each track point.

**Monitor memory.** The memory monitor in the menu bar updates every 2 seconds. If usage approaches system RAM, clear intermediate variables or spatial crop the data. The `cleanupGUI` function releases all large data arrays on close.

### SVD Cutoff Tuning Workflow (Fastest Approach)

1. Run filter once with a wide range (e.g., [1, 1000]).
2. Open the frame slider and scrub to a frame with visible bubbles.
3. Slowly increase Cutoff Start until tissue background disappears.
4. Slowly decrease Cutoff End until background noise is suppressed.
5. All adjustments after the first run are instant (cached SVD).

### Key Performance Constants

| Constant | Value | Purpose |
|----------|-------|---------|
| `DEBOUNCE_DELAY` | 50 ms | Prevents UI freezing during rapid slider movement |
| `GAUSSIAN_FIT_BATCH_SIZE` | 500 | Batch size for vectorized Gaussian fitting |
| `PROGRESS_UPDATE_INTERVAL` | 100 | Iterations between progress bar updates |
| `PLAYBACK_TIMER_PERIOD` | 50 ms | Frame advance period during playback |

---

## 17. Troubleshooting

| Symptom | Likely Cause | Solution |
|---------|-------------|----------|
| GUI won't open | Missing `.m` files or toolboxes | Ensure all files are on the MATLAB path; verify toolboxes with `ver` |
| "setDefaultParams.m is missing" warning | The file is not on the path | Non-fatal — the GUI uses safe fallback defaults. Place `setDefaultParams.m` on the path and restart. |
| "Out of Memory" during SVD | Data matrix too large for available RAM | Spatial crop first (Tab 1); close other figures (`close all`); use 64-bit MATLAB; clear workspace (`clear all`) |
| Detection count always equals Max Bubbles/Frame | Detector is saturated | Increase `Max Bubbles/Frame` until the count drops below the limit |
| SVD sliders show no change | SVD cache may be stale or frame has no signal | Click "Run Filter" to recompute; verify the displayed frame index has visible signal |
| Too many short fragmented tracks | Tracking parameters are too permissive | Increase `Min Track Length`; enable Smart Cost Matrix; reduce `Max Linking Distance`; switch to Kalman tracking |
| Localization yield < 30% of detections | FWHM mismatch or noisy filtered image | Adjust FWHM to match actual bubble PSF; reduce `Max Shift Factor`; re-tune filter cutoffs |
| Display not updating after parameter change | Processing still running or display error | Check status lamp; click frame number and press Enter to force refresh; save session and restart if needed |
| Undo button greyed out | Fewer than 2 parameter states in history | Make at least one parameter change; note that spatial crop clears undo history |
| Session fails to load | Version mismatch or missing fields | `ensureAllParamFields()` fills missing fields automatically; major structural changes may still cause errors |
| Kalman Trust Balance shows unexpected values | Panel uses theoretical gain from current parameters, not tracking data | For HKT, adjust `hk_alpha`, `hk_beta`, and `hk_v_max` in the Advanced Kalman dialog |
| "Invalid crop box" error | Crop box format incorrect | Enter as `[x y width height]` in pixels; width and height must be positive; region must be within frame bounds |
| NCC detection without PSF template | No `MB_image` loaded | Accept the auto-generate prompt to create a Gaussian PSF from the current FWHM setting |
| "No 3D numeric data found" on load | Data file doesn't contain a 3D array | Ensure the `.mat` file contains a numeric array of shape `[Nz × Nx × Nt]` |

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

## 19. Acknowledgments

This work was carried out under the supervision of **Prof. Tali Ilovitsh** at the School of Biomedical Engineering, Tel Aviv University.

The Hierarchical Kalman Tracker (`trackKalman_Advanced.m`) is based on the methodology described in:
> Taghavi, I. et al. "Ultrasound super-resolution imaging with a hierarchical Kalman tracker." *Ultrasonics*, 124, 106742, 2022. DOI: 10.1016/j.ultras.2022.106695

The SSM clutter filtering method (`SVD_SSM.m`) is based on:
> Baranger, J. et al. "Fast Thresholding of SVD Clutter Filter Using the Spatial Similarity Matrix." *IEEE TUFFC*, 70(8), 821–830, 2023. DOI: 10.1109/TUFFC.2023.3289235

The DCC-SVD method (`DCC_SVD.m`) is based on:
> Han, X. et al. "An adaptive spatiotemporal filter for ultrasound localization microscopy based on density canopy clustering." *Ultrasonics*, 144, 107446, 2024. DOI: 10.1016/j.ultras.2024.107446

The block-wise SVD method (`SVD_blockwise.m`) is based on:
> Song, P. et al. "Improved Super-Resolution Ultrasound Microvessel Imaging with Spatiotemporal Nonlocal Means Filtering and Bipartite Graph-Based Microbubble Tracking." *IEEE TBME*, 65(1), 149–167, 2017. DOI: 10.1109/TBME.2017.2703894

The benchmark in vivo rat brain dataset used for validation was provided by:
> Chavignon, A. et al. "In vivo rat brain for Ultrasound Localization Microscopy." Zenodo, 2023. DOI: 10.5281/zenodo.7883227

The radial symmetry localization algorithm (`localizeRadialSymmetry.m`) is based on:
> Parthasarathy, R. "Rapid, accurate particle tracking by calculation of radial symmetry centers." *Nature Methods*, 9(7), 724–726, 2012. DOI: 10.1038/nmeth.2071

---

*For questions, bug reports, or feature requests, please open an Issue on the GitHub repository.*
