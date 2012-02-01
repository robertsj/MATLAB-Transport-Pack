%> @file  test_inners_1d.m
%> @brief Tests the fixed source iteration in 1D
%
%       SI:  0.00000099, 1636  142.521566
% Livolant:  0.00000089,  132,  11.693956
%    GMRES:  0.00000092,   15,   6.346117
%    1.0568e+01
%    1.3411e+01
%    1.6054e+01
%    1.8535e+01
%    2.0884e+01
%    2.3122e+01
%    2.5264e+01
%    2.7323e+01
%    2.9308e+01
%    3.1227e+01
% ==============================================================================

clear

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         1);
put(input, 'inner_solver',          'Livolant');
put(input, 'inner_max_iters',       3000);
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
put(input, 'bc_left',               'vacuum');
put(input, 'bc_right',              'vacuum');


% ==============================================================================
% MATERIALS (Test one group data)
% ==============================================================================
mat = test_materials(1);

% ==============================================================================
% MESH (simple slab)
% ==============================================================================
n = 300;
xcm    = [0 1 2 3 4 5];
xfm    = n*ones(1,length(xcm)-1);
matid  = 1*ones(1,length(xcm)-1);
mesh = Mesh1D(xfm, xcm, matid);

% ==============================================================================
% FIXED SOURCE
% ==============================================================================
q_e         = Source(mesh, 1);
spectra     = [ 1.0 ];  
place       = ones(1,length(xcm)-1);
set_sources(q_e, spectra, place);


% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = GaussLegendre(64);
boundary    = BoundaryMesh(input, mesh, quadrature);
q_f         = FissionSource(state, mesh, mat);  % Not inititalized = not used.

% ==============================================================================
% SOLVE 
% ==============================================================================
solver = Fixed(input,         ...
               state,         ...
               boundary,      ...
               mesh,          ...
               mat,           ...
               quadrature,	  ...
               q_e,           ...
               q_f);
 
tic
out = solve(solver); 
toc

% ==============================================================================
% POSTPROCESS 
% ==============================================================================


% subplot(2, 1, 1)
f = flux(state, 1);
for i = 1:5
   ff(i,1) = mean(f(1+(i-1)*n : i*n)) ;
end
ff/5
% plot_flux(mesh, f)
% subplot(2, 1, 2)
% f = flux(state, 2);
% plot_flux(mesh, f)


