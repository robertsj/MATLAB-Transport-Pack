%> @file  pin_response_study.m
%> @brief Summer ANS summary -- case 2
% ==============================================================================

%clear  
warning('OFF', 'user:input')
warning('OFF', 'solver:convergence')
format longe

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-5);
put(input, 'eigen_max_iters',       100);
put(input, 'inner_tolerance',       1e-8);
put(input, 'inner_max_iters',       100);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
put(input, 'bc_left',               'reflect_r');
put(input, 'bc_right',              'reflect_r');
put(input, 'bc_top',                'reflect_r');
put(input, 'bc_bottom',             'reflect_r');

input.number_groups = 1;

mat = test_materials(1);

% Basic pin properties
p = 1.26;                 
r = [ .27 0.54];

% Materials: inside out
matid       = [2 2 1]; 

% Build the pin.
pin = PinCell(p, r, matid, 0);

% External source.
q_e         = Source(pin, 1);                  % Not initialized = not used.
spectra     = [1 0];

quadrature = CollocatedMOC(27, 1, 0, 0, 0);

%quadrature = UniformMOC(50, 1, 200);
track(pin, quadrature);


placement = [1 1 0];
set_sources(q_e, spectra, placement);

state       = State(input, pin);
q_f         = FissionSource(state, pin, mat);  % Inititalized = used.

% ==============================================================================
% SOLVE 
% ==============================================================================
ff = [ 8.925812322057839e-01
       8.594503045857023e-01     
       8.010835160443575e-01];
set_phi(state,  ff, 1);

count = 1;
phi=zeros(27, 14, 3);
for so = 0:26
for ao = 0:13
quadrature.d_order_space = so;
quadrature.d_order_angle = ao;
boundary = BoundaryTrack(input, pin, quadrature);
solver = Fixed(input,         ...
               state,         ...
               boundary,      ...
               pin,          ...
               mat,           ...
               quadrature,	...
               q_e,           ...
               q_f);
out = solve(solver); 
f = flux(state, 1); f'
phi(so+1, ao+1, 1) = f(1);
phi(so+1, ao+1, 2) = f(2);
phi(so+1, ao+1, 3) = f(3);
end
end

save('pinstudy2.mat')

return

