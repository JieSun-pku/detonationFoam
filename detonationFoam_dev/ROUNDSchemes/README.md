# libROUNDSchemes: An OpenFOAM Library of High-resolution Structure-preserving Convection Schemes

libROUNDSchemes is an open-source library for OpenFOAM. It implements ROUND schemes proposed in [1] into unstructured grids based OpenFOAM. The library can significantly reduce numerical errors with a minor increased CPU cost. Moreover, the library can offer an improved structure-preserving property that gives essentially oscillation-free solutions and preserves the structures of the passive transported scalar.

## Compilation

The libROUNDSchemes library does not require any third-party dependency.
After sourcing the OpenFOAM version, simply execute:

```
./Allwmake
```
## Usage

The library can work with any OpenFOAM solvers containing convection terms.

* The library should be linked to the solver. Add the following to your system/controlDict file:

```
libs
(
    "libROUNDSchemes.so" 
);
```

* For a general scalar field alpha, select a ROUND scheme such as ROUNDAplus in system/fvSchemes as:

```
divSchemes
{
   div(phi,alpha)      Gauss ROUNDAplus;
}
```

* For a strictly bounded scalar field such as mass fraction Y, select a ROUND scheme such as ROUNDAplus in system/fvSchemes as:

```
divSchemes
{
   div(phi,Y)      Gauss ROUNDAplus01;
}
```

* For a density-based solver, select a ROUND scheme such as ROUNDAplus in system/fvSchemes as:

```
interpolationSchemes
{
   reconstruct(rho) ROUNDAplus;
}
```

## Tutorials
The benchmark test of convection of 2D complex waves is available in _tutorial_ directory.

## Citation

If you use our library, please cite the publications describing its theory and implementation [1,2].

## Reference
- [1] Deng, X., 2023. A Unified Framework for Non-linear Reconstruction Schemes in a Compact Stencil. Part 1: Beyond Second Order. *Journal of Computational Physics*, p.112052.
- [2] Deng, X., and Xie, B., An OpenFOAM library of high-resolution structure-preserving convection schemes. (submitted)  


