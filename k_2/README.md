# k = 2: The Axisymmetric Vector Laplacian curl curl + grad div

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20%5Ctext%7Bcurl%7D%5E%7Bn*%7D_%7Brz%7D%20u%20-%20%5Ctext%7Bgrad%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%26%3D%20f%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u_%7Brz%7D%20%5Ccdot%20t%20%26%3D%200%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u_%7B%5Ctheta%7D%20%26%3D%200%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%26%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

Re-writen in its weak formulation, 

Find <img src="https://render.githubusercontent.com/render/math?math=%24(%5Csigma_h%2C%20u_h)%20%5Cin%20B_h%20%5Ctext%7Bx%7D%20C_h%24"> such that

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20(%5Csigma_h%2C%20%5Ctau_h)_r%20-%20(%5Ctext%7Bcurl%7D_%7Brz%7D%5En%20%5Ctau_h%2C%20u_h)_r%20%26%3D%200%20%5C%5C%0A%20%20%20%20%20%20%20%20(%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20%5Csigma_h%2C%20v_h)_r%20%2B%20(%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u_h%2C%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20v_h)_r%20%26%3D%20(F%2C%20v_h)_r%20%5C%5C%0A%20%20%20%20%20%20%20%20%26%20%5Cforall%20%5Ctau_h%20%5Cin%20B_h%2C%20%5Cforall%20v_h%20%5Cin%20C_h%0A%5Cend%7Baligned%7D">

## Usage

The program ...

## Usage Example
```
% add path for get_data_1() function
addpath ../data/
% define the highest mesh level
mesh = 5;
% define the nth-Fourier mode
n = 1;
% define the problem domain
pdepoly([0,1,1,0], [0,0,1,1]);
% define the equations
[u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n);
% run the program
[err_u,err_s] = weighted_HL_k_2_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,n);
```