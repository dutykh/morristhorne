<!--
  Authors: Dr. Davide Batic (Mathematics Department, Khalifa University of
  Science and Technology, Abu Dhabi, UAE) and Dr. Denys Dutykh (Mathematics
  Department, Khalifa University of Science and Technology, Abu Dhabi, UAE).
-->

# Computation of quasinormal modes for Morris-Thorne wormholes

[![Class. Quantum Grav.](https://img.shields.io/badge/Class.%20Quantum%20Grav.-41%20%282024%29%20215003-005a9c)](https://iopscience.iop.org/article/10.1088/1361-6382/ad7cb8)
[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1361--6382%2Fad7cb8-1a7f37)](https://doi.org/10.1088/1361-6382/ad7cb8)
[![arXiv](https://img.shields.io/badge/arXiv-2410.05979-b31b1b)](https://arxiv.org/abs/2410.05979)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-4c1)](LICENSE)
[![Maple](https://img.shields.io/badge/Maple-2022%2B-ff6f00)](https://www.maplesoft.com/products/Maple/)
[![MATLAB](https://img.shields.io/badge/MATLAB-Advanpix%20MCT-0076a8)](https://www.advanpix.com/)
[![QNMs Hall of Fame](https://img.shields.io/badge/QNMs-Hall%20of%20Fame-8a2be2)](https://qnms.denys-dutykh.com/)
[![GitHub](https://img.shields.io/badge/GitHub-dutykh%2Fmorristhorne-181717?logo=github)](https://github.com/dutykh/morristhorne)
[![Last commit](https://img.shields.io/github/last-commit/dutykh/morristhorne)](https://github.com/dutykh/morristhorne/commits/main)

**Authors**: Dr. Davide Batic (KU) and Dr. Denys Dutykh (KU)

**Numerical methods**: Collocation-type spectral method based on Chebyshev polynomials coupled with a polynomial eigenvalue solver

**Programming languages**: Maple (TM) and Matlab (TM)

This repository contains the computational supplemental materials for the following scientific paper:

**Reference**: D. Batic & D. Dutykh. *A unified spectral approach for quasinormal modes of Morris–Thorne wormholes*. Class. Quantum Grav., **41**, 215003, 2024

* Published version: <https://iopscience.iop.org/article/10.1088/1361-6382/ad7cb8> (DOI: `10.1088/1361-6382/ad7cb8`)
* arXiv preprint: <https://arxiv.org/abs/2410.05979>
* Code repository: <https://github.com/dutykh/morristhorne/>

**NB**: Please acknowledge this work if you use these codes.

## QNMs Hall of Fame

This work belongs to the **QNMs Hall of Fame**, a collection of studies on quasinormal modes carried out at Khalifa University and by collaborating institutions:

<https://qnms.denys-dutykh.com/>

The Hall of Fame gathers the team members, the published manuscripts, and, for each of them, the links to the journal version, the arXiv preprint and the accompanying code repository. The present paper is listed there, together with a pointer to this very repository.

## Scientific content

Morris–Thorne wormholes (also known as Bronnikov–Ellis wormholes) support damped oscillations, the quasinormal modes (QNMs), whose complex frequencies $\Omega$ are the resonances of the linearised perturbation problem. Scalar ($s = 0$), electromagnetic ($s = 1$) and gravitational ($s = 2$) perturbations all reduce to a radial wave equation with an effective potential, the spin entering through $\varepsilon = 1 - s^2$ and the angular momentum through $\lambda = L (L + 1)$.

The QNM frequencies enter the radial equation quadratically, so imposing the outgoing boundary conditions and discretising the equation produces a *quadratic* (polynomial) eigenvalue problem rather than a standard one:

```math
\left( M_0 + \mathrm{i} \Omega M_1 + \Omega^2 M_2 \right) \mathbf{a} = \mathbf{0}
```

where $\mathbf{a}$ collects the Chebyshev coefficients of the unknown eigenfunction. Three equivalent formulations of the boundary-value problem are implemented:

1. the full wormhole, both copies of the universe being mapped onto $[-1, 1]$;
2. the domain running from the throat to space-like infinity;
3. the same throat-to-infinity domain with the alternative boundary conditions introduced in the paper.

A third-order WKB computation is also provided, both as an independent check and to delimit the regime where the WKB approximation ceases to be reliable.

## Repository content

### Symbolic derivations (Maple)

| File | Content |
| --- | --- |
| `MORRIS_THORNE_WORMH_SUPPLEMENTAL_MATERIAL.mw` | Derivation of the radial ODE and of the effective potential, the change of variables encoding the QNM boundary conditions on the upper and lower copies of the universe, the analysis of the order of the zeroes of the resulting coefficients, and the verification of the limits and of the tabulated values quoted in the paper. |
| `FROM_THROAT-MORRIS-THORNE_WORMHOLE_SUPPL_MATERIAL.mw` | The same programme carried out for the formulation posed on the throat-to-infinity domain. |

### Matrix assembling (Maple)

The three worksheets below discretise the corresponding formulation with a Chebyshev collocation method and export the matrices $M_0$, $M_1$, $M_2$ of the polynomial eigenvalue problem.

| File | Formulation |
| --- | --- |
| `QNMsMorrisThorneWH-1.mw` | Full wormhole, trigonometric mapping of both universes onto $[-1, 1]$. |
| `QNMsMorrisThorneWH-2.mw` | From the throat to space-like infinity. |
| `QNMsMorrisThorneWH-3.mw` | From the throat, with the alternative boundary conditions. |

All three share the same structure and the same user-adjustable parameters, set at the top of the worksheet:

* `L` : angular momentum parameter;
* `s` : spin of the perturbation (`0` scalar, `1` electromagnetic, `2` gravitational), with `e := 1 - s^2`;
* `n` : number of Chebyshev modes, hence the size $n \times n$ of the matrices (default `200`);
* `Digits` : working precision of the symbolic-numeric assembling (default `200`).

The collocation nodes are the Chebyshev roots $x_i = \cos\left( (2i - 1) \pi / (2n) \right)$; the commented alternative in the source switches to the Chebyshev extrema of $[-1, 1]$. Each worksheet ends with three `ExportMatrix` calls writing `M0_200.mat`, `M1_200.mat`, `M2_200.mat` in ASCII Matlab format into a `data/` directory. That directory is not tracked here, since the matrices are bulky and are reproduced in a few minutes: create it (or edit the export paths) before running a worksheet.

### WKB approximation (Maple)

`QNM-WKB3.mw` implements the third-order WKB formula. It locates the peak $x_0$ of the potential $V(x) = \lambda / (1 + x^2) + \varepsilon / (1 + x^2)^2$, evaluates the derivatives of $V$ up to sixth order at that point, forms the second- and third-order corrections $\Lambda_2$, $\Lambda_3$, and returns the frequency

```math
\Omega = \sqrt{V_0 - \mathrm{i} \sqrt{-2 V_0''} \left( a + \Lambda_2 + \Lambda_3 \right)}, \qquad a = n + \frac{1}{2}
```

with `L`, `s` and the overtone number `n` set at the top of the worksheet.

### Eigenvalue solution (Matlab)

| File | Content |
| --- | --- |
| `findeigsdp.m` | Double-precision solution. Reads the three exported matrices, calls `polyeig(M0, 1i*M1, M2)` and plots the computed spectrum in the complex $\Omega$ plane. |
| `findeigshp.m` | Multiple-precision counterpart. Sets `mp.Digits(200)`, reads the matrices with `mp.read` and performs the same polynomial eigenvalue computation in extended precision. |

`findeigsdp.m` looks for the matrices in `data/`, whereas `findeigshp.m` reads them from the current folder; adjust either the export paths or the reading paths so that the two agree.

The multiple-precision version relies on the [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/); adjust the `addpath` line at the top of `findeigshp.m` to point to your local installation. Collocation discretisations of QNM problems are notoriously ill-conditioned, which is the reason for the 200 working digits used throughout: the double-precision spectrum is best regarded as indicative, and the digits to be trusted are those confirmed by the extended-precision run.

## Typical workflow

1. Open one of the `QNMsMorrisThorneWH-*.mw` worksheets in Maple, choose `L`, `s` and `n`, and execute it. The three matrices are written to `data/`.
2. Run `findeigsdp.m` in Matlab for a quick double-precision look at the spectrum.
3. Run `findeigshp.m` for the multiple-precision computation, which delivers the digits quoted in the paper.
4. Optionally, run `QNM-WKB3.mw` for the third-order WKB estimate of the same frequency and compare.

## Requirements

* Maple 2022 or later (the worksheets are saved in that format; they use `LinearAlgebra` and `PDEtools`).
* Matlab with the built-in `polyeig`.
* Advanpix Multiprecision Computing Toolbox, for `findeigshp.m` only.

## Citation

```bibtex
@article{Batic2024,
  author  = {Batic, Davide and Dutykh, Denys},
  title   = {A unified spectral approach for quasinormal modes of {Morris--Thorne} wormholes},
  journal = {Classical and Quantum Gravity},
  volume  = {41},
  number  = {21},
  pages   = {215003},
  year    = {2024},
  doi     = {10.1088/1361-6382/ad7cb8}
}
```

## License

The codes are distributed under the terms of the GNU Lesser General Public License, version 2.1. See `LICENSE`.
