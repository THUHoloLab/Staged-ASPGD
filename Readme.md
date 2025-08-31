# ASPGD Flattop Beam Shaping  (MATLAB)

This repository contains a MATLAB script that simulates **amplitude- and phase-aberrated** beam shaping using a spatial phase (SLM-like) mask and an **Adaptive Stochastic Parallel Gradient Descent (ASPGD)** optimizer. The goal is to produce a **flattop (super-Gaussian-like)** intensity at the Fourier plane while compensating Zernike-mode phase errors and an elliptical amplitude distortion.

---

## Contents

* `main.m` — the full simulation script (the code you pasted)
* *(expected helper functions; see “Required helper functions” below)*

  * `Fx_Zernike.m`
  * `Fx_Guassian.m` *(Gaussian field generator; note the spelling in the script)*
  * `Fx_STD.m` *(PSF sharpness metric)*
  * `Fx_MSD.m` *(mean-squared difference metric)*
  * `Fx_Structure.m` *(structural similarity/shape metric)*

---

## What the script does

1. **Defines a Fourier optics model** (scalar, paraxial) with a circular pupil:

   * Units are **millimeters** throughout (wavelength given in mm).
   * Builds input and output sampling grids consistent with an f-plane FFT model.

2. **Constructs aberrations**

   * **Phase**: linear combination of the first `Zern_AP` Zernike modes with user-defined coefficients.
   * **Amplitude**: separable elliptical Gaussian envelope.

3. **Computes an “ideal” flattop phase**

   * Uses a stationary-phase-like separable phase `φ_in(x)+φ_in(y)` whose Fourier magnitude approximates a square/flat profile (via error-function primitives).

4. **Evaluates pre-optimization performance**

   * Shows input intensity/phase and the unoptimized output intensity.

5. **Runs ASPGD optimization in two stages**

   * **Stage 1 (PSF sharpening)**: drives the PSF toward a delta-like spot using a sharpness metric (`Fx_STD`).
   * **Stage 2 (target shaping)**: minimizes an **α-scaled MSD** to a flattop target plus a **structure** term, using Zernike coefficients as control variables.

6. **Plots results**

   * Target vs. achieved output, before/after line profiles, and cost histories.

---

## Requirements

* **MATLAB** R2021a+ (earlier versions likely fine)
* **Toolboxes**

  * *Statistics and Machine Learning Toolbox* (for `binornd`)
* **Helper functions** (provide these in the MATLAB path):

  * `Fx_Zernike(R)`: returns a Zernike stack `Zern_N(:,:,k)` over the pupil of radius `R`.
  * `Fx_Guassian(R0)`: returns the circularly symmetric Gaussian field amplitude with 1/e² radius `R0` on the **input grid**.
  * `Fx_STD(I)`: scalar sharpness metric (e.g., std/variance or Tenengrad; you choose).
  * `Fx_MSE(I, T)`: mean squared error between image `I` and target `T`.
  * `Fx_Structure(I, T)`: structure/shape similarity in \[0,1] (1 = identical structure).



---

## How to run

1. Place `main.m` and all helper `.m` files in the same folder (or add them to the MATLAB path).
2. Open MATLAB and run:

   ```matlab
   clear; clc; close all;
   main
   ```
3. Inspect the figures:

   * **Input intensity/phase**
   * **Ideal phase & output**
   * **Before/after optimization** images and **horizontal line profiles**
   * **Cost curves** for each stage

---

## Key parameters (edit at the top of the script)

| Name             | Meaning                                              | Default                           |
| ---------------- | ---------------------------------------------------- | --------------------------------- |
| `lambda`         | Wavelength (mm)                                      | `1064e-6` (1064 nm)               |
| `NA`             | Numerical aperture                                   | `0.45`                            |
| `f`              | Focal length (mm)                                    | `10`                              |
| `dx0`            | Input-plane sampling (mm/pixel)                      | `12.5e-3`                         |
| `R_0`            | 1/e² **intensity** radius of the Gaussian input (mm) | `1.5`                             |
| `L_flat`         | Flattop half-width in the output plane (mm)          | `0.1`                             |
| `Zern_AP`        | Number of phase aberration modes injected            | `12`                              |
| `Zern_Level`     | Number of optimized Zernike modes                    | `65`                              |
| `Iteration_N`    | Iterations (per stage)                               | `200` (stage 1), `2000` (stage 2) |
| `Noiseamplitude` | SPGD probe amplitude                            | `0.1` (stage 1), `0.3` (stage 2)  |
| `gama`           | Structure-term weight                                | `100`                             |

> **Grids and sizes**
>
> * Pupil radius: `R = NA * f`
> * Input FOV: `L0 = 2*R` and `N = ceil(L0/dx0)`
> * Output sampling: `dx1 = lambda * f / L0` (FFT scaling)

---

## Algorithm details (ASPGD / SPSA-style)

* At each iteration, draw **Rademacher** perturbations `Δc ∈ {±Noiseamplitude}^K` over Zernike coefficients.
* Evaluate **Add/Sub** costs with `c ± Δc`.
* Form the stochastic gradient estimate

  $$
  g \approx \frac{J(c+\Delta c)-J(c-\Delta c)}{2\,\|\Delta c\|_\infty^2}\,\Delta c
  $$
* Apply **Adam-like** momentum/variance normalization:
  `momentum ← β₁·momentum + (1−β₁)·g`
  `veloc    ← β₂·veloc    + (1−β₂)·(g.^2)`
  `Δc_adam  ← momentum ./ sqrt(veloc + ε)`
* Update coefficients: `c ← c + η · Δc_adam` (η = `StepLength`, possibly scheduled).

> ⚠️ **Important**: In the posted script `Zern_Level` is same in both stage 1 and 2 for simplicity.

**Stage-specific costs**

* **Stage 1**: `J = Fx_STD(PSF)` (maximize sharpness ⇒ minimize negative sharpness or similar).
* **Stage 2**:

`J = α · MSD(I, I_t) + γ · (1 − Structure(I,I_t))^0.1`
  where `α = 1/Fx_MSD(F3d, I_t)` at stage start (normalization).

---

## Figures produced

* **Ideal design**: stationary phase `φ_in(x,y)`, corresponding output, input field, and a line profile.
* **Before/after**: target vs. optimized output, pre-optimization output, profiles through the center row.
* **Cost histories**: `CostV` (stage 1) and `CostVM` (stage 2).


---

## How to adapt

* **Change the target**: Edit construction of `I_t` to any pattern (e.g., ring, multi-spot). And calculate the corresponding phase `phi_in`
* **Add amplitude control**: Present code varies only **phase** via Zernike modes. To include amplitude masks, extend the parameterization and the forward model.
* **Different metrics**: Swap `Fx_STD`, `Fx_MSD`, and `Fx_Structure` for your preferred choices (e.g., SSIM, MTF-based metrics).

---

## Minimal checklist before running

* [ ] Provide all helper `.m` files listed above.
* [ ] Confirm `binornd` is available (Statistics Toolbox) or replace with `sign(randn(size))`.
* [ ] Verify units (mm) for your `f`, `NA`, `dx0`, `R_0`, and desired `L_flat`.

---

## Citation / attribution

Please cite the following paper if the code is used in your research work. Thank you!

Churan Han, Liangcai Cao, Dun Liu, Hao Tu, and Qiaofeng Tan, "High uniformity flattop beam shape correction with complex amplitude aberration of the incidence," Optics and Lasers in Engineering 195, 109275 (2025).

https://doi.org/10.1016/j.optlaseng.2025.109275

Contact: clc@tsinghua.edu.cn

---



