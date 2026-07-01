# DG_Clarke

MATLAB implementation of the **DG-Clarke** method for derivative-free nonsmooth optimization using discrete gradients and MatCSG search directions.

This repository contains the main solver, a driver for the finite-max test problems, the objective evaluator, the test-environment file, and the hitlist files used for benchmarking.

The related repository

```text
https://github.com/GS1400/DG-Clarke-TEminmax
```

contains the supplementary numerical material for the paper, including `suppMat/suppMat.pdf`, the construction of the `TE.mat` test environment, and additional scripts for checking and summarizing the finite-max problems.

## Contents

The repository is organized as follows:

```text
DG_Clarke/
    DG_Clarke.m
    driver_DG_Clarke.m
    README.md
    LICENSE
    TE/
        TE.mat
        TE_realML_minmax_eval.m
    HIT/
        RMLMMAX001.mat
        RMLMMAX002.mat
        ...
        RMLMMAX080.mat
```

## Main files

### `DG_Clarke.m`

The file

```text
DG_Clarke.m
```

contains the implementation of the DG-Clarke method. The solver uses discrete-gradient approximations and an inner DG-MatCSG procedure.

The two direction modes are:

```matlab
tune.dir = 1;   % steepest discrete-gradient direction
tune.dir = 2;   % MatCSG direction
```

The main call is:

```matlab
[x,outAll] = DG_Clarke(fun,x0,tune,ST);
```

where `fun` is the objective-function handle, `x0` is the starting point, `tune` contains algorithmic parameters, and `ST` contains stopping and printing parameters.

### `driver_DG_Clarke.m`

The file

```text
driver_DG_Clarke.m
```

runs the solver on one problem from the finite-max test environment and compares the two direction modes `dir = 1` and `dir = 2`.

For example:

```matlab
results = driver_DG_Clarke(1);
```

runs the method on problem `RMLMMAX001`.

The driver loads the corresponding problem from `TE/TE.mat`, loads the hitlist file from `HIT/`, runs both direction modes, and produces a comparison plot of `log10(q_f)` versus the number of function evaluations.

## Test environment

The folder

```text
TE/
```

contains the MATLAB file

```text
TE/TE.mat
```

which stores the structure `TE`. This structure contains the 80 real-data and synthetic finite-max nonsmooth test problems.

In MATLAB:

```matlab
load('TE/TE.mat','TE')
TE
```

The structure `TE` contains general information such as:

```matlab
TE.description
TE.nproblem
TE.ndataset
TE.nloss
TE.created
TE.seed
TE.problem
TE.names
```

The field `TE.problem` contains the individual problem structures:

```matlab
TE.problem.RMLMMAX001
TE.problem.RMLMMAX002
...
TE.problem.RMLMMAX080
```

For example:

```matlab
prob = TE.problem.RMLMMAX001;
```

returns the data for problem `RMLMMAX001`, including the problem name, dimension, data set, loss type, number of finite-max terms, bounds, starting points, and objective-function handle.

The objective evaluator for the real-data finite-max problems is:

```text
TE/TE_realML_minmax_eval.m
```

## Starting point

The benchmark driver uses the stored initial point

```matlab
TE.problem.(pname).points.xr
```

for all reported runs.

The test environment also contains other reference points, such as `points.x0` and `points.x1`, but the driver uses `points.xr` to match the hitlist generation and benchmarking setup.

## Hitlist files

The folder

```text
HIT/
```

contains the hitlist files for all 80 finite-max problems:

```text
HIT/RMLMMAX001.mat
HIT/RMLMMAX002.mat
...
HIT/RMLMMAX080.mat
```

Each file contains a MATLAB structure `hitlist`.

For example:

```matlab
load('HIT/RMLMMAX034.mat','hitlist')
hitlist
```

returns a structure such as:

```matlab
hitlist = 

  struct with fields:

        xopt: [8×1 double]
        fopt: 0.3347
    normxopt: 1.3370
```

Here, `xopt` is the best available reference point, `fopt` is the corresponding best available reference value, and `normxopt` is the norm of the reference point.

The value `fopt` is used as the benchmark target value in the stopping ratio

```matlab
q_f = (f_best - fopt)/(f0 - fopt).
```

This target value is used only for benchmarking and is not assumed to be available in practical optimization runs.

## Git LFS

The file

```text
TE/TE.mat
```

is large and is tracked using Git LFS.

After cloning the repository, make sure Git LFS is installed:

```bash
git lfs install
git lfs pull
```

To check that the large file is tracked correctly:

```bash
git lfs ls-files
```

The output should include:

```text
TE/TE.mat
```

## Basic MATLAB usage

After cloning the repository, add all folders to the MATLAB path:

```matlab
addpath(genpath(pwd));
```

Load the test environment:

```matlab
load('TE/TE.mat','TE')
```

Access one problem:

```matlab
prob = TE.problem.RMLMMAX001;
```

Evaluate the objective function at a point `x`:

```matlab
f = prob.funf(x);
```

Load the hitlist for the same problem:

```matlab
load('HIT/RMLMMAX001.mat','hitlist')
fopt = hitlist.fopt;
xopt = hitlist.xopt;
```

Run the driver on one problem:

```matlab
results = driver_DG_Clarke(1);
```

Run the driver on all 80 problems:

```matlab
for p = 1:80
    fprintf('\nRunning problem %03d\n',p);
    results = driver_DG_Clarke(p);
end
```

## Output of the driver

For each problem, the driver runs:

```matlab
tune.dir = 1;
tune.dir = 2;
```

and reports the final objective value, the final benchmark ratio, the number of logged function evaluations, and the stopping flag.

The driver also saves comparison plots and result files. Generated files such as plots, summaries, and result files are ignored by Git through `.gitignore`.

## Related repository

The companion repository

```text
https://github.com/GS1400/DG-Clarke-TEminmax
```

## Funding

This project was funded by the **Austrian Science Fund (FWF)** under the project

**Derivative-free methods for nonsmooth optimization (DFNO)**

**Grant DOI:** [10.55776/PAT2747625](https://doi.org/10.55776/PAT2747625)


contains the supplementary numerical material for the paper

**An Approximate Conjugate Subgradient Algorithm with Matrix Parameter for Derivative-Free Nonsmooth Optimization Problems**

by **Morteza Kimiaei**, **Saman Babaie--Kafaki**, and **Zohre Aminifard**.

That repository includes:

```text
suppMat/suppMat.pdf
TEminmax/TE.mat
hitlist/RMLMMAX001.mat
...
hitlist/RMLMMAX080.mat
```

It also describes how the test environment was generated, how the finite-max test problems are organized, and how the benchmark target values were obtained.

## Citation

If you use this code, the finite-max test environment, `TE.mat`, or the hitlist files, please cite the accompanying paper:

M. Kimiaei, S. Babaie-Kafaki, and Z. Aminifard,
*An Approximate Conjugate Subgradient Algorithm with Matrix Parameter for Derivative-Free Nonsmooth Optimization Problems*,
arXiv:2606.29084, 2026.
https://arxiv.org/abs/2606.29084

```bibtex
@misc{kimiaei_babaiekafaki_aminifard_2026_dgclarke,
  title         = {An Approximate Conjugate Subgradient Algorithm with Matrix Parameter for Derivative-Free Nonsmooth Optimization Problems},
  author        = {Kimiaei, Morteza and Babaie-Kafaki, Saman and Aminifard, Zohre},
  year          = {2026},
  eprint        = {2606.29084},
  archivePrefix = {arXiv},
  primaryClass  = {math.OC},
  url           = {https://arxiv.org/abs/2606.29084},
  note          = {DG-Clarke MATLAB implementation and supplementary numerical material}
}
```

## License

See the file

```text
LICENSE
```

for licensing information.
