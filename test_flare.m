%> @file  test_flare.m
%> @brief 2-d flare model.
%clear classes
% 
% map = [3 2 2 2 3 2 2 1 0
%        0 2 2 2 2 2 2 1 0 
%        0 2 2 2 2 2 1 1 0 
%        0 2 2 2 2 2 1 0 0 
%        0 2 2 2 3 1 1 0 0 
%        0 2 2 2 1 1 0 0 0 
%        0 2 1 1 1 0 0 0 0 
%        0 1 1 0 0 0 0 0 0 
%        0 0 0 0 0 0 0 0 0];
   
   
 map = [2, 0, 1, 0, 1, 0, 1, 2,-1 
        -1,1, 0, 1, 0, 1, 0, 2,-1 
        -1,0, 1, 0, 1, 0, 1, 2,-1 
        -1,1, 0, 1, 0, 1, 2, 2,-1 
        -1,0, 1, 0, 1, 0, 2,-1,-1
        -1,1, 0, 1, 0, 1, 2,-1,-1  
        -1,0, 1, 2, 2, 2,-1,-1,-1 
        -1,2, 2, 2,-1,-1,-1,-1,-1
       -1,-1,-1,-1,-1,-1,-1,-1,-1]
map = map+1;
%   double precision :: temp_d1(4)  = (/1.4493e+00, 1.4479e+00, 1.4494e+00, 1.3200e+00/)
%   double precision :: temp_d2(4)  = (/3.8070e-01, 3.7080e-01, 3.6760e-01, 2.7720e-01/)
%   double precision :: temp_r1(4)  = (/2.5000e-02, 2.5800e-02, 2.6200e-02, 0.02576220/)
%   double precision :: temp_a2(4)  = (/1.0420e-01, 1.2000e-01, 1.1910e-01, 7.1596e-02/)
%   double precision :: temp_f1(4)  = (/7.9000e-03, 6.9000e-03, 6.0000e-03, 0.0000e+00/)
%   double precision :: temp_f2(4)  = (/1.6920e-01, 1.7450e-01, 1.6250e-01, 0.0000e+00/)
%   double precision :: temp_s12(4) = (/1.5100e-02, 1.4800e-02, 1.4700e-02, 2.3106e-02/)
% 0.2300    0.2302    0.2300    0.2525 s1
% 0.8756    0.8990    0.9068    1.2025 s2
mat = Materials(2, 2);

% fresh
set_diff_coef(mat, 1, 1,    1.4493e+00);
set_diff_coef(mat, 1, 2,    3.8070e-01);
set_sigma_t(mat, 1, 1,      0.2300);
set_sigma_t(mat, 1, 2,      0.8756);
set_sigma_s(mat, 1, 1, 1,   0.2300 - 2.5000e-02);
set_sigma_s(mat, 1, 2, 1,   1.5100e-02);
set_sigma_s(mat, 1, 2, 2,   0.8756 - 1.0420e-01);
set_nu_sigma_f(mat, 1, 1,   7.9000e-03);
set_nu_sigma_f(mat, 1, 2,   1.6920e-01);
set_chi(mat, 1, 1,          1.0);

% once
set_diff_coef(mat, 2, 1,    1.4479e+00);
set_diff_coef(mat, 2, 2,    3.7080e-01);
set_sigma_t(mat, 2, 1,      0.2302);
set_sigma_t(mat, 2, 2,      0.8990);
set_sigma_s(mat, 2, 1, 1,   0.2302 - 2.5800e-02);
set_sigma_s(mat, 2, 2, 1,   1.4800e-02);
set_sigma_s(mat, 2, 2, 2,   0.8990 - 1.2000e-01);
set_nu_sigma_f(mat, 2, 1,   7.9000e-03);
set_nu_sigma_f(mat, 2, 2,   1.6920e-01);
set_chi(mat, 2, 1,          1.0);

% twice
set_diff_coef(mat, 3, 1,    1.4479e+00);
set_diff_coef(mat, 3, 2,    3.7080e-01);
set_sigma_t(mat, 3, 1,      0.2300);
set_sigma_t(mat, 3, 2,      0.9068);
set_sigma_s(mat, 3, 1, 1,   0.2300 - 2.6200e-02);
set_sigma_s(mat, 3, 2, 1,   1.4700e-02);
set_sigma_s(mat, 3, 2, 2,   0.9068 - 1.1910e-01);
set_nu_sigma_f(mat, 3, 1,   6.0000e-03);
set_nu_sigma_f(mat, 3, 2,   1.6250e-01);
set_chi(mat, 3, 1,          1.0);

input = Input();
put(input, 'node_width',    21);
put(input, 'mixing_factor', 0.9);
put(input, 'albedo_single', 0.2);
put(input, 'albedo_double', 0.6); 

solver = FLARE(input, mat, map);


[s,k]=solve(solver);
plot_power(solver);