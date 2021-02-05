# k = 0: The Neumann Problem for the Axisymmetric Poisson Equation

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20-%20%5Ctext%7Bdiv%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20u%20%26%20%3D%20f%20%26%26%5Ctext%7B%20in%20%7D%20%5COmega%2C%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20u%20%5Ccdot%20n%20%26%20%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

Re-writen in its weak formulation, 

Find <img src="https://render.githubusercontent.com/render/math?math=%24Q_h%20%5Cin%20A_h%24"> such that

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20(%5Ctext%7Bgrad%7D%5E%5Ctext%7Bn%7D_%7Brz%7DQ_h%20u%2C%20%5Ctext%7Bgrad%7D%5E%5Ctext%7Bn%7D_%7Brz%7D%20v_h)_%7BL%5E2_r(%5COmega)%7D%20%26%3D%20(%5Ctext%7Bgrad%7D%5E%5Ctext%7Bn%7D_%7Brz%7Du%2C%20%5Ctext%7Bgrad%7D%5E%5Ctext%7Bn%7D_%7Brz%7D%20v_h)_%7BL%5E2_r(%5COmega)%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%26%20%5Cforall%20v_h%20%5Cin%20A_h%20%26%0A%5Cend%7Baligned%7D">

Note: <img src="https://render.githubusercontent.com/render/math?math=%24A_h%20%3D%20%5C%7B%20u%20%5Cin%20(H_r(%5Ctext%7Bgrad%7D%5En%2C%20%5COmega)%5C%7D%24"> for a given mesh <img src="https://render.githubusercontent.com/render/math?math=%24%5COmega%24">

## Usage


### Syntax
```
[err] = weighted_HL_k_0_p1_e(u,grad_u_r,grad_u_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z
```

```
[err,grad_err,max_err] = weighted_HL_k_0_e(u,grad_u_r,grad_u_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z)
```

### Inputs
`f` - given function  
`grad_f_r` - gradient(f) with respect to r  
`grad_f_z` - gradient(f) with respect to z  
`gd,sf,ns` - outputs of pdepoly specifying domain  
`mesh` - max mesh level  
`n` - n-th Fourier mode  
`u` - exact solution function  
`grad_u_r` - gradient(u) with respect to r  
`grad_u_z` - gradient(u) with respect to z  

### Outputs
`err` - array of L2 errors for mesh levels corresponding to indices  
`grad_err` - array of L2 gradient errors for mesh levels corresponding to indices  
`max_err` - array of max errors for mesh levels corresponding to indicies  

## Example
```
% add path for get_data_7() function
addpath ../data/
% define the highest mesh level
mesh = 5;
% define the nth-Fourier mode
n = 1;
% define the problem domain
pdepoly([0,1,1,0], [0,0,1,1]);
% define the equations
[f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data7();
% run the program
% use the first for P1 (first order) or the second for P2 (second order)
[err] = weighted_HL_k_0_p1_e(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z);
[err,grad_err,max_err] = weighted_HL_k_0_p2_e(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z);
```

Return to [main](../README.md) page
