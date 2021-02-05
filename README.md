# Finite Element Methods to Approximate Hodge Laplacian Problems on an Axisymmetric Domain

## Table of Contents

* [Introduction](#introduction)
* [Equations](#equations)
* [Definitions](#definitions)
* [Usage](#usage)
* [Requirements](#requirements)

## Introduction

This repository contains efficient finite element methods to approximate Hodge Laplacian problems on an axisymmetric domain.

An axisymmetric problem is one defined on a three-dimensional domain that is symmetric with respect to an axis. By using cylindrical coordinates and a Fourier series decomposition with respect to <img src="https://render.githubusercontent.com/render/math?math=\theta"> , the three-dimensional problem can be reduced to a sequence of two-dimensional problems, which are easier to solve and significantly reduces computation time. The two-dimensional problems are posed in weighted function spaces with weight r.

The finite element method (FEM) is a numerical technique for approximating solutions to complex differential equations. A large portion of differential equations cannot be solved using analytical techqniques; rather thay must be approximated. The FEM is an ideal candidate for approximating solutions due to its well developed theory, adaptability, and accuracy.

## Equations
<!-- https://jsfiddle.net/8ndx694g/ Converts LaTex equations to rendered URLs -->
The finite element methods approximate the solution to the following weighted mixed formulation of the abstract Hodge Laplacian.

Find <img src="https://render.githubusercontent.com/render/math?math=(\sigma , u) \in V^{k-1} x V^k"> such that:

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20(%5Csigma%20%2C%20%5Ctau)_%7BL%5E2_r(%5COmega)%7D%20-%20(d%5E%7Bk-1%7D%20%5Ctau%20%2C%20u)_%7BL%5E2_r(%5COmega)%7D%20%26%3D%200%2C%20%26%26%20%5Ctext%7B%20for%20all%20%7D%20%5Ctau%20%5Cin%20V%5E%7Bk-1%7D%2C%20%5C%5C%0A%20%20%20%20(d%5E%7Bk-1%7D%5Csigma%20%2C%20v)_%7BL%5E2_r(%5COmega)%7D%20%2B%20(d%5Ek%20u%2C%20d%5Ek%20v)_%7BL%5E2_r(%5COmega)%7D%20%26%3D%20(f%2Cv)_%7BL%5E2_r(%5COmega)%7D%20%2C%20%26%26%20%5Ctext%7B%20for%20all%20%7D%20v%20%5Cin%20V%5Ek.%0A%5Cend%7Baligned%7D">

With <img src="https://render.githubusercontent.com/render/math?math=%24k%20%3D%200%2C1%2C2%2C3%24">.

See the [next](#definitions) section for more details.

Thus, the four equations corresponding with <img src="https://render.githubusercontent.com/render/math?math=%24k%20%3D%200%2C1%2C2%2C3%24"> are as follows:

### k = 0: The Neumann Problem for the Axisymmetric Poisson Equation

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20-%20%5Ctext%7Bdiv%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20u%20%26%20%3D%20f%20%26%26%5Ctext%7B%20in%20%7D%20%5COmega%2C%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20u%20%5Ccdot%20n%20%26%20%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

See [here](k_0/README.md) for more

### k = 1: The Axisymmetric Vector Laplacian curl curl + grad div

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20-%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20%5Ctext%7Bdiv%7D%5E%7Bn*%7D_%7Brz%7D%20%2B%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bcurl%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20u%20%26%3D%20f%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20(%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20u)_%7Brz%7D%20%5Ccdot%20t%20%26%3D%200%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20(%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20u)_%7B%5Ctheta%7D%20%26%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u_%7Brz%7D%20%5Ccdot%20n%20%26%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

See [here](k_1/README.md) for more

### k = 2: The Axisymmetric Vector Laplacian curl curl + grad div

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20%5Ctext%7Bcurl%7D%5E%7Bn*%7D_%7Brz%7D%20u%20-%20%5Ctext%7Bgrad%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%26%3D%20f%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u_%7Brz%7D%20%5Ccdot%20t%20%26%3D%200%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u_%7B%5Ctheta%7D%20%26%3D%200%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%26%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

See [here](k_2/README.md) for more

### k = 3: The Dirichlet Problem for the Axisymmetric Poisson Equation

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20-%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20%5Ctext%7Bgrad%7D%5E%7Bn*%7D_%7Brz%7D%20u%20%26%3D%20f%20%26%26%5Ctext%7B%20in%20%7D%20%5COmega%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u%20%26%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

See [here](k_3/README.md) for more


### Other Equations

#### Weighted Lowest Order Nedelec Space

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%5Coverrightarrow%7B%5Ctext%7Bcurl%7D%7D_%7Brz%7D%20%5Ctext%7Bcurl%7D_%7Brz%7D%20u%20%20%2B%20u%20%26%20%3D%20F%20%26%26%5Ctext%7B%20in%20%7D%20%5COmega%2C%20%5C%5C%0A%20%20%20%20%5Ctext%7Bcurl%7D_%7Brz%7D%20u%20%26%20%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%5Cend%7Baligned%7D">

#### Weighted Fourier Lowest Order Nedelec and P1 Space (Bh)

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%5Ctext%7Bcurl%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20u%20%20%2B%20u%20%26%20%3D%20F%20%26%26%5Ctext%7B%20in%20%7D%20%5COmega%2C%20%5C%5C%0A%20%20%20%20(%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20u)_%7Brz%7D%20%5Ccdot%20t%20%20%26%20%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1%2C%5C%5C%0A%20%20%20%20(%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20u)_%7B%5Ctheta%7D%20%20%26%20%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%5Cend%7Baligned%7D">

#### Lowest Order Raviart Thomas Space

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20-%20%5Ctext%7Bgrad%7D%20(%5Ctext%7Bdiv%7D%20u)%20%2B%20u%20%26%3D%20F%2C%5C%5C%0A%20%20%20%20%20%20%20%20%5Ctext%7Bdiv%7D%20u%20%26%3D%200%20%26%26%5Ctext%7Bon%20%7D%20%5Cpartial%5COmega.%0A%20%20%20%20%5Cend%7Baligned%7D">

#### Weighted Fourier Lowest Order Raviart Thomas Space (Ch)

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20-%20%5Ctext%7Bgrad%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%2B%20u%20%26%3D%20F%2C%5C%5C%0A%20%20%20%20%20%20%20%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%26%3D%200%20%26%26%5Ctext%7Bon%20%7D%20%5Cpartial%5COmega.%0A%5Cend%7Baligned%7D%0A">

## Definitions

<!-- d^k-->
We let <img src="https://render.githubusercontent.com/render/math?math=d^k"> be defined in the following way:

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0A%20%20%20%20d%5E0%20v%20%26%3D%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20v%2C%5C%5C%0A%20%20%20%20d%5E1%20v%20%26%3D%20%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20v%2C%5C%5C%0A%20%20%20%20d%5E2%20v%20%26%3D%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20v%2C%5C%5C%0A%20%20%20%20d%5E3%20v%20%26%3D%200%2C%0A%5Cend%7Balign*%7D">

<!-- V^k-->
and let <img src="https://render.githubusercontent.com/render/math?math=V^k"> be the Hilbert space associated with each,

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign*%7D%0A%20%20%20%20V%5E0%20%26%3D%20H_r(%5Ctext%7Bgrad%7D%5En%2C%20%5COmega)%2C%5C%5C%0A%20%20%20%20V%5E1%20%26%3D%20H_r(%5Ctext%7Bcurl%7D%5En%2C%20%5COmega)%2C%5C%5C%0A%20%20%20%20V%5E2%20%26%3D%20H_r(%5Ctext%7Bdiv%7D%5En%2C%20%5COmega)%2C%5C%5C%0A%20%20%20%20V%5E3%20%26%3D%20L%5E2_r(%5COmega).%0A%5Cend%7Balign*%7D">

Furthermore, the grad, curl, and div formulas for the n-th Fourier mode are as follows:
<!-- grad^n_rz-->
<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%5Ctext%7Bgrad%7D%5E%7Bn%7D_%7Brz%7D%20v%20%26%3D%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cpartial_r%20v%20%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20-%20%5Cfrac%7Bn%7D%7Br%7D%20v%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cpartial_z%20v%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%2C%5C%5C%0A%20%20%20%20%20%20%20%20%5Ctext%7Bcurl%7D%5E%7Bn%7D_%7Brz%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_r%20%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_%7B%5Ctheta%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_z%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%20%26%3D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20-(%20%5Cfrac%7Bn%7D%7Br%7D%20v_z%20-%20%5Cpartial_z%20v_%7B%5Ctheta%7D%20)%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cpartial_z%20v_r%20-%20%5Cpartial_r%20v_z%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7Bn%20v_r%20%2B%20v_%7B%5Ctheta%7D%7D%7Br%7D%20%2B%20%5Cpartial_r%20v_%7B%5Ctheta%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%2C%5C%5C%0A%20%20%20%20%20%20%20%20%5Ctext%7Bdiv%7D%5E%7Bn%7D_%7Brz%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_r%20%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_%7B%5Ctheta%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_z%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%20%26%3D%20%5Cpartial_r%20v_r%20%2B%20%5Cfrac%7Bv_r%20-%20n%20v_%7B%5Ctheta%7D%7D%7Br%7D%20%2B%20%5Cpartial_z%20v_z%20.%0A%20%20%20%20%5Cend%7Baligned%7D">

<!-- grad^n*_rz-->

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%5Ctext%7Bgrad%7D%5E%7Bn*%7D_%7Brz%7D%20v%20%26%3D%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cpartial_r%20v%20%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7Bn%7D%7Br%7D%20v%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cpartial_z%20v%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%2C%5C%5C%0A%20%20%20%20%20%20%20%20%5Ctext%7Bcurl%7D%5E%7Bn*%7D_%7Brz%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_r%20%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_%7B%5Ctheta%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_z%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%20%26%3D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7Bn%7D%7Br%7D%20v_z%20-%20%5Cpartial_z%20v_%7B%5Ctheta%7D%20%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cpartial_z%20v_r%20-%20%5Cpartial_r%20v_z%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7B-n%20v_r%20%2B%20v_%7B%5Ctheta%7D%7D%7Br%7D%20%2B%20%5Cpartial_r%20v_%7B%5Ctheta%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%2C%5C%5C%0A%20%20%20%20%20%20%20%20%5Ctext%7Bdiv%7D%5E%7Bn*%7D_%7Brz%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%5B%20%7B%5Cbegin%7Barray%7D%7Bcc%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_r%20%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_%7B%5Ctheta%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20v_z%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Barray%7D%20%7D%20%5Cright%5D%20%26%3D%20%5Cpartial_r%20v_r%20%2B%20%5Cfrac%7Bv_r%20%2B%20n%20v_%7B%5Ctheta%7D%7D%7Br%7D%20%2B%20%5Cpartial_z%20v_z%20.%0A%20%20%20%20%5Cend%7Baligned%7D">

TODO: For more information, see thesis

## Requirements

* MATLAB

## James Madison University Honors Capstone Project
Author: Nicole Stock

Faculty Research Advisor: Dr. Minah Oh
