%> @file  test_fixed_1d.m
%> @brief Tests the fixed source iteration in 1D
%
%> Finish me.
% ==============================================================================

clear

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'inner_solver',          'Livolant');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
input.number_groups = 2;


% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);


% ==============================================================================
% MESH (The simple slab reactors from Mosher's thesis)
% ==============================================================================

% Assembly coarse mesh edges
base   = [ 1.1580 4.4790 7.8000 11.1210 14.4420 15.6000 ]; 
% and fine mesh counts.
basef  = [ 1 2 2 2 2 1 ]*4; 
% Several such assemblies to make the total coarse mesh definition
xcm    = [ 0.0  base  base+15.6  base+15.6*2 base+15.6*3 ...
           base+15.6*4 base+15.6*5 base+15.6*6 ];
xfm    = [ basef basef basef basef basef basef basef ];

% Assembly types
assem_A    = [ 1 2 3 3 2 1 ];
assem_B    = [ 1 2 2 2 2 1 ];
assem_C    = [ 1 2 4 4 2 1 ];
assem_D    = [ 1 4 4 4 4 1 ];

% Cores 1, 2 and 3
core_1 = [ assem_A assem_B assem_A assem_B assem_A assem_B assem_A ];
core_2 = [ assem_A assem_C assem_A assem_C assem_A assem_C assem_A ];
core_3 = [ assem_A assem_D assem_A assem_D assem_A assem_D assem_A ];

mesh = Mesh1D(xfm, xcm, core_1);

% ==============================================================================
% FIXED SOURCE
% ==============================================================================
q_e         = Source(mesh, 2);
spectra     = [ 0.0 1.0
                0.0 0.0 ];  
bases  = [ 1 2 2 2 2 1 ];     
place    = [ bases bases bases bases bases bases bases ];
set_sources(q_e, spectra, place);


% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = GaussLegendre(8);
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
% f = flux(state, 1);
% plot_flux(mesh, f)
% subplot(2, 1, 2)
% f = flux(state, 2);
% plot_flux(mesh, f)


