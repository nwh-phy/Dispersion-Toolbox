# q-EELS Dispersion Fitting Toolbox

A MATLAB-based interactive toolbox for **visualizing and fitting dispersion relations** from momentum-resolved Electron Energy Loss Spectroscopy (q-EELS) data.

Features include Lorentz peak fitting on individual spectra, automated batch fitting across all momentum channels, dispersion curve extraction, and publication-quality figure export.

---

## Quick Start

```matlab
% Launch with preprocessed eq3D data
interactive_qe_browser("path/to/eq3D.mat");

% Launch with previously saved operation history
interactive_qe_browser("path/to/eq3D.mat", "path/to/op_history.mat");
```

**Input format**: `eq3D.mat` file containing a 3D intensity matrix (`eq3D`) and momentum resolution (`dq`).

---

## Features

- **Four-panel interactive GUI**: q-E heatmap, single spectrum viewer, comparison map, dispersion plot
- **Lorentz peak fitting**: Drude-Lorentz model with power-law background on individual spectra
- **Automated batch fitting**: fit all q-channels at once with initial peak guesses
- **Dispersion extraction**: manual peak picking, branch splitting, and model fitting
- **Linewidth analysis**: Γ(q) and Γ(E) visualization in a separate window
- **Signal processing**: Gaussian/Savitzky-Golay smoothing, Wiener denoising, Lucy-Richardson deconvolution
- **On-axis separation**: symmetric off-axis referencing to isolate excitation signals from ZLP tails
- **Background subtraction**: Power-law, ExpPoly3, and Pearson models
- **Export**: PNG (300 dpi) + vector PDF with selectable aspect ratios
- **Operation history**: full undo/redo, save/load session history

---

## Architecture

```
┌──────────────────────────────────────────────────┐
│            interactive_qe_browser.m              │
│        (GUI: visualization · fitting · export)   │
├──────────┬───────────────┬───────────────────────┤
│ Data I/O │  Analysis     │  Fitting              │
├──────────┼───────────────┼───────────────────────┤
│load_qe_  │build_         │fit_loss_function      │
│dataset   │comparison_qe  │  (Lorentz + BG)       │
│make_qe_  │               │fit_quasi2d_plasmon    │
│struct    │               │  (dispersion model)   │
└──────────┴───────────────┴───────────────────────┘
```

---

## Core Functions

| Function | Description |
|----------|-------------|
| `load_qe_dataset` | Load `eq3D.mat` and build standardized q-E data structure with automatic q-axis calibration |
| `make_qe_struct` | Construct q-E struct from intensity matrix, energy axis, and dq |
| `build_comparison_qe` | On-axis signal separation using symmetric off-axis reference bands |
| `fit_loss_function` | Single-spectrum Lorentz fitting with power-law background model |
| `fit_quasi2d_plasmon` | Dispersion relation fitting with user-defined analytical model |

---

## GUI Layout

| Panel | Position | Description |
|-------|----------|-------------|
| **q-E Heatmap** | Top-left | Momentum-energy map; click to select q-channel |
| **Spectrum** | Top-right | Single spectrum with Lorentz fit overlay |
| **Comparison** | Bottom-left | Normalized off-axis signal map |
| **Dispersion** | Bottom-right | Extracted dispersion points + fitted curve |

---

## Workflow

### 1. Load & Configure
```matlab
interactive_qe_browser("eq3D.mat");
% Set q range, energy range, smoothing parameters
```

### 2. Preprocess
```
Optional: enable Denoise / BG Sub / Deconv / Area Norm
Set reference q-range → click Build Views for on-axis separation
```

### 3. Fit Spectra
```
Click on q-E map to select a channel
→ Enter peak guesses (e.g., "800,2500") → Fit Spectrum
→ Accept Fit to add peaks to dispersion plot
Or: Auto Fit ω(q) for batch fitting across all channels
```

### 4. Fit Dispersion
```
Split branches → Fit Model → extract model parameters
Show Γ → analyze linewidth trends
```

### 5. Export
```
Select aspect ratio → 📷 Export → choose panels
→ PNG (300 dpi) + PDF (vector) output
```

---

## Fitting Models

### Loss Function (Single Spectrum)

$$S(E) = B_0 \cdot E^{-\alpha} + \sum_{i=1}^{N} \frac{A_i \cdot E \cdot \Gamma_i}{(E^2 - \omega_{p,i}^2)^2 + E^2 \cdot \Gamma_i^2}$$

- Power-law background models the ZLP tail
- Each Lorentz peak is parameterized by ($\omega_p$, $\Gamma$, $A$)

### Dispersion Model

The default dispersion model follows [da Jornada et al., *Nat. Commun.* **11**, 1013 (2020)]:

$$E(q) = \sqrt{\frac{A \cdot |q|}{\varepsilon_{\text{bg}} + \rho_0 \cdot |q|}}$$

This can be replaced with any user-defined analytical dispersion relation.

---

## Requirements

- **MATLAB R2024b+**
- Optimization Toolbox
- Signal Processing Toolbox

---

## File Structure

```
├── interactive_qe_browser.m   # Main interactive GUI
├── load_qe_dataset.m          # Data loading
├── make_qe_struct.m           # q-E struct constructor
├── build_comparison_qe.m      # On-axis signal separation
├── fit_loss_function.m        # Lorentz + power-law fitting
├── fit_quasi2d_plasmon.m      # Dispersion model fitting
└── references/                # Reference literature
```

---

## License

This project is provided for academic and research use.

## Citation

If you use this toolbox in your research, please cite:

> da Jornada, F. H., Xian, L., Rubio, A., & Louie, S. G. (2020). Universal slow plasmons and giant field enhancement in atomically thin quasi-two-dimensional metals. *Nature Communications*, **11**, 1013.
