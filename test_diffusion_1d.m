%> @file  test_diffusion_1d.m
%> @brief Tests 1D diffusion stuff.
%
%> Finish me.
% ==============================================================================

clear

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'dimension',         1);
put(input, 'number_groups',     1);
put(input, 'inner_solver',      'SI');
put(input, 'inner_print_out',   1);
put(input, 'outer_print_out',   1);

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(1);

% ==============================================================================
% MESH (The simple slab reactors from Mosher's thesis)
% ==============================================================================


% Several such assemblies to make the total coarse mesh definition
xcm    = [ 0.0  10];
xfm    = [   10 ];
mt     = [  1 ];


mesh = Mesh1D(xfm, xcm, mt);

% ==============================================================================
% FIXED SOURCE
% ==============================================================================
q_e         = Source(mesh, 1);
spectra     = [ 1.0 ];     
place       = [ 1 ];
set_sources(q_e, spectra, place);


% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = GaussLegendre(8);
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

f = flux(state, 1); 

diffop = DiffusionOperator(input, mat, mesh);
M = diffop.get_1g_operator(1);


