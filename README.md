## Introduction

This code is an implementation of the subproblem solver that is presented in the article

1. ```Marvin Severitt and Paul Manns. Efficient solution of discrete subproblems arising in integer optimal control with total variation regularization. INFORMS Journal on Computing, 35 (4), 2023.``` 

for subproblems that appear in the sequential linear integer programming algorithm (SLIP) for TV-regularized integer
optimal control problems, which is presented in the article 

2. ```Sven Leyffer and Paul Manns. Sequential linear integer programming for integer optimal control with total variation regularization. ESAIM: Control, Optimisation and Calculus of Variations (ESAIM: COCV) 28 (66), 2022.```

Currently, the implementation only supports uniform discretizations of a one-dimensional control domain ```(a,b)```.


## Targets for CMake

* ```trs4slip```: Main library
* ```trs4slip_test```: Unit tests
* ```trs4slip_python_sdist```: Source package for python bindings
* ```trs4slip_python_whl```: Compiled (wheel) package for python bindings

## Dependencies

For the python bindings, you will need a python interpreter with the `build`, `numpy`, and `dev` packages being installed. Currently, python>=3.10 is required. You can also change the version requirement in the CMakeLists.txt to 3.8.

## Install the C++ library (paths relative to project root)

```console
$ mkdir -p build
$ cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
$ sudo make install # <-- system-wide install
```

## Test if the algorithm works correctly after building the C++ library (paths relative to project root)

```console
$ cd build
$ make test
```

## Install python package (paths relative to project root)

(load your anaconda workspace, etc.)

```console
$ mkdir -p build
$ cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make trs4slip_python_whl
$ cd ../python/dist
$ ls
```
Now, there is a file `SOME_COMPLICATED_FILENAME.whl` in this folder. The python module can be installed now

```console
$ pip install SOME_COMPLICATED_FILENAME.whl
```

## Run python example (paths relative to project root)

```console
$ cd examples
$ python run_slip_algorithm.py
```

## Usage of the C++ and python interface

The python package provides a function `run` that encapsulates
the C++-function `trs4slip_astar`, which hints at the fact that a variant
of the A\* algorithm is used to solve the subproblems. 

The python package provides a function `run_top` that encapsulates
the C++-function `trs4slip_top`, which hints at the fact that a variant
of the topological sorting-based algorithm is used to solve the subproblems. 

The C++ interface is as follows.

```cpp
void trs4slip_astar(
    int32_t * x_next_out,
    double const *c,
    int32_t const *x,
    int32_t const *bangs,
    int const delta,
    int const N,
    int const M,
    double * vert_costs_buffer,
    int32_t * vert_layer_buffer,
    int32_t * vert_value_buffer,
    int32_t * vert_prev_buffer,
    int32_t * vert_remcap_buffer,
    const int boundcon,
    double lbound,
    double rbound  
);

```
* ```x_next_out```: used to store (and return) the solution computed by the subproblem solver. The user needs to provide an integer array of size N if N denotes the number of discretization intervals for the control. The quantity corresponds to ```x + d``` in [1].
* ```c```: used to store the cost coefficients (see (TR-IP) on p. 4 in [1]). The user needs to provide the cost as a double array of size N. The quantity corresponds to ```c``` in [1].
* ```x```: used to store the current iterate. The user needs to provide the current iterate as an integer array of size N if N denotes the number of discretization intervals for the control. The integers need to be elements of bangs (below). The quantity corresponds to ```x``` in [1].
* ```bangs```: used to store the feasible control realizations. The user needs to provide an integer array of size M if N denotes the number of feasible control realizations. The quantity corresponds to &Xi; in [1].
* ```delta```: used to store the trust-region radius as an integer. It is the L1-norm-value scaled by N.
* ```N```: number of discretization intervals as integer.
* ```M```: number of discrete realizations as integer.
* ```vert_costs_buffer```: buffer for internal purposes in order to serve frequent reallocation of memory. The user needs to provide a double array of size ```N*M*(Delta_max + 1) + 2```, where ```Delta_max``` is the maximum trust-region radius delta that is ever used by the caller.
* ```vert_layer_buffer```: buffer for internal purposes in order to serve frequent reallocation of memory. The user needs to provide an integer array of size ```N*M*(Delta_max + 1) + 2```, where ```Delta_max``` is the maximum trust-region radius delta that is ever used by the caller.
* ```vert_value_buffer```: buffer for internal purposes in order to serve frequent reallocation of memory. The user needs to provide an integer array of size ```N*M*(Delta_max + 1) + 2```, where ```Delta_max``` is the maximum trust-region radius delta that is ever used by the caller.
* ```vert_prev_buffer```: buffer for internal purposes in order to serve frequent reallocation of memory. The user needs to provide an integer array of size ```N*M*(Delta_max + 1) + 2```, where ```Delta_max``` is the maximum trust-region radius delta that is ever used by the caller.
* ```vert_remcap_buffer```: buffer for internal purposes in order to serve frequent reallocation of memory. The user needs to provide an integer array of size ```N*M*(Delta_max + 1) + 2```, where ```Delta_max``` is the maximum trust-region radius delta that is ever used by the caller.
* ```boundcon```: flag that indicates if Dirichlet boundary conditions are supplied for x (implying that jumps to these boundary values are also taken into account in the objective value).
* ```lbound```: left boundary value (only read if boundcon not false).
* ```rbound```: right boundary value (only read if boundcon not false).

```cpp
void trs4slip_top(
    int32_t * x_next_out,
    double const *c,
    int32_t const *x,
    int32_t const *bangs,
    double const * switchingcost,
    int const delta,
    int const N,
    int const M,
    int actbounds,
    double leftbound,
    double rightbound,
    double switchingcostleft,
    double switchingcostright
);

```
* ```x_next_out```: used to store (and return) the solution computed by the subproblem solver. The user needs to provide an integer array of size N if N denotes the number of discretization intervals for the control. The quantity corresponds to ```x + d``` in [1].
* ```c```: used to store the cost coefficients (see (TR-IP) on p. 4 in [1]). The user needs to provide the cost as a double array of size N. The quantity corresponds to ```c``` in [1].
* ```x```: used to store the current iterate. The user needs to provide the current iterate as an integer array of size N if N denotes the number of discretization intervals for the control. The integers need to be elements of bangs (below). The quantity corresponds to ```x``` in [1].
* ```bangs```: used to store the feasible control realizations. The user needs to provide an integer array of size M if N denotes the number of feasible control realizations. The quantity corresponds to &Xi; in [1].
* ```switchingcost```: used to store switching cost coefficients for the jump points between the interior intervals. The user needs to provide double array of size N - 1. Setting all values to 1 gives the total variation.
* ```delta```: used to store the trust-region radius as an integer. It is the L1-norm-value scaled by N.
* ```N```: number of discretization intervals as integer.
* ```M```: number of discrete realizations as integer.
* ```actbounds```: flag that indicates if Dirichlet boundary conditions are supplied for x (implying that jumps to these boundary values are also taken into account in the objective value).
* ```leftbound```: left boundary value (only read if actbounds not false).
* ```rightbound```: right boundary value (only read if actbounds not false).
* ```switchingcostleft```: switching cost coefficient for the jump point to the left boundary value (only read if actbounds not false).
* ```switchingcostright```: switching cost coefficient for the jump point to the right boundary value (only read if actbounds not false).

The python interface is as follows.

```python
def run(np.int32_t[::1] x_next not None,
        double[::1] c not None,
        np.int32_t[::1] x not None,
        np.int32_t[::1] bangs not None,
        int delta,
        double[::1] vert_costs_buffer not None,
        np.int32_t[::1] vert_layer_buffer not None,
        np.int32_t[::1] vert_value_buffer not None,
        np.int32_t[::1] vert_prev_buffer not None,
        np.int32_t[::1] vert_remcap_buffer not None,
        int boundcon,
        double lbound,
        double rbound
        ): 
```

```python
def run_top(np.int32_t[::1] x_next not None, # output argument
        double[::1] c not None,
        np.int32_t[::1] x not None,
        np.int32_t[::1] bangs not None,
        double[::1] switchingcost not None,
        int delta,
        int actbounds,
        double leftbound,
        double rightbound,
        double switchcostleft,
        double switchcostright
        ):
```

The variables are the same as in the C++ case but N and M are deduced from the given numpy arrays.
In python, the allocation of the working memory arrays may look as follows.

```python
vert_costs_buffer = np.empty(N*M*(Delta_max + 1) + 2)
vert_layer_buffer = np.empty(N*M*(Delta_max + 1) + 2, dtype=np.int32)
vert_value_buffer = np.empty(N*M*(Delta_max + 1) + 2, dtype=np.int32)
vert_prev_buffer = np.empty(N*M*(Delta_max + 1) + 2, dtype=np.int32)
vert_remcap_buffer = np.empty(N*M*(Delta_max + 1) + 2, dtype=np.int32)
```
