%> @file  test_pin_cell.m
%> @brief Tests the eigenvalue solver on a 2D pin cell problem.
%  4.905678061750672e-01
%  4.908480133215041e-01 < ref sn
%  Reference: No pin   -- 4.909795825291901e-01 (scatter,
%  6.508113024900996e-01; 6.563226932826230e-01
%             with pin -- 5.513339209625725e-01     4.086689341553977e-01
% sc:  7.429774497100994e-01     5.363325792754957e-01
%> Finish me. 
% ==============================================================================

clear  
warning('OFF', 'user:input')
warning('OFF', 'solver:convergence')

sn = 0;

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-5);
put(input, 'eigen_max_iters',       100);
put(input, 'inner_tolerance',       0.0001);
put(input, 'inner_max_iters',       100);
put(input, 'inner_solver',          'SI');
put(input, 'livolant_free_iters',   3);
put(input, 'livolant_accel_iters',  3);
put(input, 'bc_left',               'reflect');
%put(input, 'bc_right',              'reflect');
%put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'reflect');

input.number_groups = 1;

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(1);


% ==============================================================================
% MESH (A simple pin-cell)
% ==============================================================================

% Basic pin properties
p = 1.26;                 
r = [ .27 0.54];

% Materials: inside out
matid       = [2 2 1]; 

% Build the pin.
pin = PinCell(p, r, matid, 0);

% External source.
q_e         = Source(pin, 1);                  % Not initialized = not used.
spectra     = [1 1];      

if sn == 1
    % Mesh the pin.
    meshify(pin, 30);
    placement = ones(number_cells(pin), 1);
    % plot
    figure(1)
    plot_mesh(pin);
    map = reshape(mesh_map(pin, 'REGION'), number_cells(pin), 1);
    m1 = map==1;
    m2 = map==2;
    m3 = map==3;
    quadrature = LevelSymmetric(8);
    boundary   = BoundaryMesh(input, pin, quadrature);
    set_sources_mesh(q_e, spectra, placement);     
else
    quadrature = CollocatedMOC(9, 1, 0);
    
    %quadrature = UniformMOC(50, 1, 200);
    track(pin, quadrature);
    clf
    %Plotting tracks
     figure(2)
    plot_tracks(pin, 1);
   %verify_tracks(pin)
    boundary = BoundaryTrack(input, pin, quadrature);
    placement = [1 1 1]; 
    set_sources(q_e, spectra, placement);    
end

% 
% ==============================================================================
% SETUP 
% ==============================================================================

phi = [8.907616871619879e-01     8.579010883562891e-01     8.051670051698306e-01]';
state       = State(input, pin);
%set_phi(state, phi, 1);

q_f         = FissionSource(state, pin, mat);  % Inititalized = used.
%initialize(q_f);

% build_reflect_index(boundary);
% build_side_index(boundary);
% Sweep 
% source  = [1; 1];
% equation = SCMOC(pin, mat);
% initialize(boundary, 1);
% sweeper = SweepMOC(input, pin, mat, quadrature, boundary, equation);
% phi = sweep(sweeper, source, 1);


% ==============================================================================
% SOLVE 
% ==============================================================================
solver = Fixed(input,         ...
                     state,         ...
                     boundary,      ...
                     pin,          ...
                     mat,           ...
                     quadrature,	...
                     q_e,           ...
                     q_f);
 
tic
    out = solve(solver); 
toc

%phiref = [4.909795825291901e-01];
% phiref = [5.513339209625725e-01     4.086689341553977e-01]; %4.909795825291901e-01
% 
 f = flux(state, 1); f=f';
if sn == 1
    if pin.d_number_regions == 2
    phi = [ mean(f(m1)) mean(f(m2)) ]
    elseif pin.d_number_regions == 3
    phi = [ mean(f(m1)) mean(f(m2)) mean(f(m3)) ]
    else
    phi = mean(f)
    end
else
    phi = f
end


return

% VTK.savevtk(state, 2, pin, 'pin.vtk')
% % Get the flux
% f = flux(state, 1);
% % and plot it
% subplot(2,1,1)
% plot_flux(pin, f)
% axis square
% % Get the flux
% f = flux(state, 2);
% % and plot it
% subplot(2,1,2)
% plot_flux(pin, f)
% axis square