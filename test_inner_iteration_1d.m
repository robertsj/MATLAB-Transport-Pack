%> @file  test_inners_1d.m
%> @brief Tests the fixed source iteration in 1D
%
%       SI: 684
% Livolant:  66
%    GMRES:
% ==============================================================================

clear

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         1);
put(input, 'inner_solver',          'Livolant');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);

% ==============================================================================
% MATERIALS (Test one group data)
% ==============================================================================
mat = test_materials(1);

% ==============================================================================
% MESH (simple slab)
% ==============================================================================

xcm    = 0:0.625:100;
xfm    = 4*ones(1,length(xcm)-1);
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
quadrature  = GaussLegendre(4);
boundary    = Boundary(input, mesh, quadrature);
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
f(1:10)
% plot_flux(mesh, f)
% subplot(2, 1, 2)
% f = flux(state, 2);
% plot_flux(mesh, f)


