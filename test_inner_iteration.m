%> @file  test_inner_iteration.m
%> @brief Tests the inner iteration.
%
%> For the defaults, should converge in 29 iterations at a rate of 0.69361896.
%> Min flux = 1.6285, Max flux = 4.5219
% ==============================================================================

clear

% Get the default input.
input               = Input();
input.tol_inner     = 1e-6;
input.solver_type   = 'Livolant';

put(input, 'livolant_free_iters',  3);
put(input, 'livolant_accel_iters', 3);

% One material, one group, c = 0.9
mat         = test_materials(1);

% Simple uniformly meshed square
mesh        = test_mesh(1);

quadrature  = LevelSymmetric(8);
% State and boundary
state       = State(input, mesh);
boundary    = Boundary(input, mesh, quadrature);

% Uniform source.
q_e         = Source(mesh, 1);
spectra     = [ 1.0 ];      
placement   = [ 1 ];    
set_sources(q_e, spectra, placement);

% Empty fission source.
q_f = FissionSource(mesh, mat);

% Make the inner iteration.
if strcmp(input.solver_type, 'SourceIteration')
    inner = SourceIteration();              
elseif strcmp(input.solver_type, 'Livolant')
    inner = Livolant();
else
    
end
setup(  inner,  ...
        input,   	...
        state,    	...
        boundary,   ...
        mesh,     	...
        mat,       	...
        quadrature, ...
        q_e,     	...
        q_f);
    

% Solve the problem
tic
[err, it] = solve(inner, 1); 
t = toc;
disp(['  norm residual = ', num2str(err)])
disp(['           time = ', num2str(t)])
disp(['         sweeps = ', num2str(it)])
disp(['time per sweep  = ', num2str(t/it)])
% Get the flux
f = flux(state, 1);
% and plot it
plot_flux(mesh, f);


% Flux should be (quarter map)
%     1.1372    1.6066    1.7923    1.9010    1.9450
%     1.6066    2.3072    2.5879    2.7520    2.8200 
%     1.7923    2.5879    2.9344    3.1327    3.2203
%     1.9010    2.7520    3.1327    3.3470    3.4427 
%     1.9450    2.8200    3.2203    3.4427    3.5552    
