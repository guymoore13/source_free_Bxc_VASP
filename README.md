# source_free_Bxc_VASP

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
    - `LPMCURRENT`: set to `True` to write $j_p$ to the `JPARAMAG` file.

## Visualization

## How to cite

## License
