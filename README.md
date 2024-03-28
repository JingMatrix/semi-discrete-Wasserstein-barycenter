## Design goal

In the Euclidean plane,
calculate Wasserstein barycenter for one
absolutely continuous marginal measure and other discrete marginal measures.

### Current state
1. The absolutely continuous measure could be the uniform measure on any (convex) polygon in the Euclidean plane.
2. For semi-discrete optimal transport, which coincides with the calculation of Wasserstein barycenter of an absolutely
continuous measure and a discrete measure, this code gives quick result for the discrete measure of size 1200.
In the author's laptop, it takes around 15s.

### Todo lists
- [ ] Compute for absolutely continuous measure that is not uniform. There is no essential difference, we just need to add integration code.
- [ ] Support non-convex polygon, either to improve the power diagram generation code, either we fall back to the first task.
- [ ] Implement a special linear programming routine for current problem, to reduce the memory usage.

### Difficulties

1. The semi-discrete solvers doesn't respect the geometric meaning of a partition. It sometimes push a vertex out of the support, and this will in turn fail the numerical solver itself.

## Compile and use

### dependencies

- cmake
- CGAL
- gsl
- glpk

### build

```sh
mkdir build
cmake -B build
make -C build -j 4
```
### cmake issue of CGAL on Arch Linux

On Arch Linux, the cmake config file [Installation/lib/cmake/CGAL/CGALConfig.cmake](https://github.com/CGAL/cgal/blob/master/Installation/lib/cmake/CGAL/CGALConfig.cmake) might fail to set CGAL_ROOT due to symbolic linking problem.

One may thus need to change its lines
```cmake
set(CGAL_ROOT ${CGAL_CONFIG_DIR})
get_filename_component(CGAL_ROOT "${CGAL_ROOT}" DIRECTORY)
```

to

```cmake
set(CGAL_ROOT ${CGAL_CONFIG_DIR})
get_filename_component(CGAL_ROOT "${CGAL_ROOT}" REALPATH)
```

### use

To use the code, see `src/test` as an example.
The above build command will build an executable `build/test`, for which one can do some experiments.
Either the user can provide his own data, or he may use the `python3` script `data/randomdata.py` to generate
some random data.
Feel free to modify this script for his own interest.

## Copyrights

All rights and permissions are reserved.
