%> @file  pin_response_study.m
%> @brief Tests the eigenvalue solver on a 2D pin cell problem.
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
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');

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
ff = [ 1.433381242517740e+00     
      1.467147540150487e+00    
       1.544447166145315e+00];
set_phi(state,  ff, 1);

count = 1;
phi=zeros(27, 14, 3);
% for so = 0:26
% for ao = 0:13
% quadrature.d_order_space = so;
% quadrature.d_order_angle = ao;
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
% end
% end



return

