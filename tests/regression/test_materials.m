%> @file  test_materials.m
%> @brief Return two group, four material benchmark data.
%
%> Reference: 
%>   Ilas and Rahnema, "A Heterogeneous Coarse Mesh Transport
%>     Method", Transport Theory and Statistical Physics, 32, 5-7, 
%>     pp445-472 (2003)
function mat = test_materials(tester)

% Create the Materials object.
mat = Materials(2, ...  % groups
                4);     % materials

% ---------------------------
% Material 1: Water
% ---------------------------

% Total
mat = set_sigma_t(mat, 1, 1, 0.1890);       % (obj, matid, g, value);
mat = set_sigma_t(mat, 1, 2, 1.4633);

% Fission
mat = set_nu_sigma_f(mat, 1, 1, 0.0);       % Note, default is zero
mat = set_nu_sigma_f(mat, 1, 2, 0.0);
mat = set_chi(mat, 1, 1, 0.0);
mat = set_chi(mat, 1, 2, 0.0);

% Scattering
mat = set_sigma_s(mat, 1, 1, 1, 0.1507);    % 1 <- 1
mat = set_sigma_s(mat, 1, 1, 2, 0.0000);    % 1 <- 2
mat = set_sigma_s(mat, 1, 2, 1, 0.0380);    % 2 <- 1
mat = set_sigma_s(mat, 1, 2, 2, 1.4536);    % 2 <- 2

% ---------------------------
% Material 2: Fuel I
% ---------------------------

% Total
mat = set_sigma_t(mat, 2, 1, 0.2263);       % (obj, matid, g, value);
mat = set_sigma_t(mat, 2, 2, 1.0119);

% Fission
mat = set_nu_sigma_f(mat, 2, 1, 0.0067);
mat = set_nu_sigma_f(mat, 2, 2, 0.1241);
mat = set_chi(mat, 2, 1, 1.0);
mat = set_chi(mat, 2, 2, 0.0);

% Scattering
mat = set_sigma_s(mat, 2, 1, 1, 0.2006);    % 1 <- 1
mat = set_sigma_s(mat, 2, 1, 2, 0.0000);    % 1 <- 2
mat = set_sigma_s(mat, 2, 2, 1, 0.0161);    % 2 <- 1
mat = set_sigma_s(mat, 2, 2, 2, 0.9355);    % 2 <- 2

% ---------------------------
% Material 3: Fuel II
% ---------------------------

% Total
mat = set_sigma_t(mat, 3, 1, 0.2252);       % (obj, matid, g, value);
mat = set_sigma_t(mat, 3, 2, 0.9915);

% Fission
mat = set_nu_sigma_f(mat, 3, 1, 0.0078);
mat = set_nu_sigma_f(mat, 3, 2, 0.1542);
mat = set_chi(mat, 3, 1, 1.0);
mat = set_chi(mat, 3, 2, 0.0);

% Scattering
mat = set_sigma_s(mat, 3, 1, 1, 0.1995);    % 1 <- 1
mat = set_sigma_s(mat, 3, 1, 2, 0.0000);    % 1 <- 1
mat = set_sigma_s(mat, 3, 2, 1, 0.0156);    % 1 <- 1
mat = set_sigma_s(mat, 3, 2, 2, 0.9014);    % 1 <- 1

% ---------------------------
% Material 4: Fuel II + Gd
% ---------------------------

% Total
mat = set_sigma_t(mat, 4, 1, 0.2173);       % (obj, matid, g, value);
mat = set_sigma_t(mat, 4, 2, 1.0606);

% Fission
mat = set_nu_sigma_f(mat, 4, 1, 0.0056);
mat = set_nu_sigma_f(mat, 4, 2, 0.0187);
mat = set_chi(mat, 4, 1, 1.0);
mat = set_chi(mat, 4, 2, 0.0);

% Scattering
mat = set_sigma_s(mat, 4, 1, 1, 0.1902);	% 1 <- 1
mat = set_sigma_s(mat, 4, 1, 2, 0.0000);    % 1 <- 1
mat = set_sigma_s(mat, 4, 2, 1, 0.0136);    % 1 <- 1
mat = set_sigma_s(mat, 4, 2, 2, 0.5733);    % 1 <- 1

% ---------------------------
% FINALIZE
% ---------------------------

% This sets the scattering bounds, which can eliminate a few operations.
mat = finalize(mat);

% ---------------------------
% VERIFY
% ---------------------------

tests = {'mat.sigma_t(1, 1) == 0.1890', ...
         'mat.nu_sigma_f(2, 2) == 0.1241', ...
         'mat.sigma_s(3, 2, 1) == 0.0156', ...
         'number_materials(mat) == 4'};
tester.run_tests(tests);


end