# `source_free_Bxc_VASP`

## Introduction

This code repository contains the necessary patch files for the implementation of the source-free constraint on $B_{xc}$ ( $\nabla \cdot B_{xc} = 0$ ) in VASP, which we explore in [this paper](preprint_url). Our implementation follows from the work of [Sharma et al.](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01049), with the key distinction that we do not include the magnetization rescaling (i.e. $s = 0$ ).

Our implementation involves the use of parallel three-dimensional Fast Fourier Transforms in order to solve the Poisson equation for the sources in the original $B_{xc}$. Therefore it is quite efficient and requires little additional computational cost compared to noncollinear VASP.

## Installation

To apply the relevant code patches, simply run the following `bash` script in your terminal:

```
$ ./apply_patch.sh
```

The installation script will walk you through the relevant installation steps, including:

1) The choice of VASP version to which you will be applying the patch. The current implementations include:
  * VASP version 6.2.1
  * VASP version 5.4.4

2) The flavor of patch: "minimal" or "full":
  * The <ins>__minimal__</ins> patch requires fewer patches to the original source-code, and simply applies the constraint directly, such that if one runs the resulting `vasp_ncl`, the resulting $B_{xc}$ will be free of monopoles ( $\nabla \cdot B_{xc} = 0$ ).
  * The <ins>__full__</ins> patch provides the following VASP `INCAR` flags to toggle the source-free constraint, as well as the I/O of various fields relevant to the source-free correction
    - `LSOURCEFREE`: toggle the $\nabla \cdot B_{xc} = 0$ constraint. `LSOURCEFREE=False` will result in the original `vasp_ncl` behavior.
    - `LVXC`: set to `True` to write $v_{xc}$ and $B_{xc}$ to the `XCPOT` file.
    - `LSOURCEPOT`: set to `True` to write $\nabla \phi$ to the `SOURCEPOT` file.
    - `LAXC`: set to `True` to write $A_{xc}$ to the `AXCPOT` file.
    - `LPMCURRENT`: set to `True` to write $j_p$ to the `JPARAMAG` file. Based on error messages from VASP, the paramagnetic current can only be calculated when `KPAR=1`.

To keep with the conventions of VASP, most of the volumetric data files (`XCPOT`, `SOURCEPOT`, and `AXCPOT`) are written in the format of the noncollinear [CHGCAR](https://www.vasp.at/wiki/index.php/CHGCAR#Noncollinear_magnetism) ( $\rho$, $m_x$, $m_y$, $m_z$ ). For the `XCPOT`, the first, ${\rho}$-component is the self-consistent $\nabla \cdot B_{xc}$, and the latter three indices contain the $x$, $y$, and $z$ components of $B_{xc}$.

## Visualization

A Jupyter notebook has been included, `visualization/Bxc_visualize.ipynb`. This minimal notebook provides the user with the means to plot $B_{xc}$, and other fields output by the patch using `plotly` and `pymatgen`, as well as other Python packages.

As an example, below is a visualization of the magnetization (from `CHGCAR`), $m(r)$, as a vector field, as well as $B_{xc}(r)$ (from `XCPOT`) as streamlines (`plotly`'s `Streamtube`), for the source-free ground-state of Mn<sub>3</sub>ZnN.

<img src="visualization/example_files/Mn3ZnN/Mn3ZnN_source_free.png" alt="fields" width="600"/>

## How to cite

## License

This repository is licensed under and MIT License shown in `LICENSE`.
