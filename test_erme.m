%> @file  test_erme.m
%> @brief Demo of ERME and solvers.

% ==============================================================================
% Server (here, a DB)
% ==============================================================================
clear
% db_input = Input();
% put(db_input, 'rf_db_name', 'two_group_3x3.h5');
% %put(db_input, 'rf_db_name', 'homog.h5');
% server = ResponseDB(db_input);
% read_response(server);

% on the fly
path(path,'./test');
inputrf = Input();
put(inputrf, 'number_groups',         2);
put(inputrf, 'inner_tolerance',       1e-9);
put(inputrf, 'inner_max_iters',       100);
put(inputrf, 'outer_tolerance',       1e-9);
put(inputrf, 'outer_max_iters',       10);
put(inputrf, 'print_out',             0);
put(inputrf, 'quad_order',            2);
so=1; 
put(inputrf, 'rf_max_order_space',       so);
put(inputrf, 'rf_max_order_azimuth',      3);
put(inputrf, 'rf_max_order_polar',        0);
put(inputrf, 'rf_order_space',       so);
put(inputrf, 'rf_order_azimuth',      3);
put(inputrf, 'rf_order_polar',        0);
put(inputrf, 'rf_number_nodes', 3);
put(inputrf, 'owner',           'rf');

% Materials
mat = mat_2g_4m();
% Pins
[pin1, pin2, pin3, pin4] = test_pins();
pins = {pin1, pin2, pin3, pin4};
% Assemblies
pin_map1 = [1 1 1 ; 1 1 1 ; 1 1 1]; % 3x3 UO2 pins
assem1   = Assembly(pins, pin_map1); meshify(assem1);
pin_map2 = [1 1 1 ; 1 2 1 ; 1 1 1]; % 3x3 UO2 with Gd type I central
assem2   = Assembly(pins, pin_map2); meshify(assem2);
pin_map3 = [1 1 1 ; 1 3 1 ; 1 1 1]; % 3x3 UO2 with Gd type II central
assem3   = Assembly(pins, pin_map3); meshify(assem3);
pin_map4 = [4 4 4 ; 4 4 4 ; 4 4 4]; % Homogeneous Gd type I fuel
assem4   = Assembly(pins, pin_map4); meshify(assem4);
assemblies = {assem1, assem2, assem3};
server = ResponseServer(inputrf,mat,assemblies);
%[R, F, A, L] = get_responses(server, 1);


% ==============================================================================
% ERME setup
% ==============================================================================
input = Input();
put(input, 'bc_left',           'reflect');
put(input, 'bc_right',          'reflect');
put(input, 'bc_bottom',         'reflect');
put(input, 'bc_top',            'reflect');
put(input, 'inner_tolerance',       1e-10);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-10);
put(input, 'outer_max_iters',       20);
put(input, 'rf_order_space',   so); % Eventually, these will be the
put(input, 'rf_order_azimuth',  3); % orders extracted from a DB
put(input, 'rf_order_polar',    0); 
put(input, 'number_groups',     2);
put(input, 'dimension',         2);
put(input, 'owner',           'erme');

elements = [3 2 1 2 1
                 2 1 2 1 2
                 1 2 1 2 1
                 2 1 2 1 2
                 1 2 1 2 1];     
problem = ERME(input, server, elements);
solver = ERME_Picard(input, problem);
tic, solve(solver); toc
% 1.319782 0.891295 1.297046 | 1.278287 <-- 1x1 reflective
% 2g, assem1, vac 1 side, k = 0.551020, 3x3 assem1 k =  0.781512
% 3x3 assem1 + assem2 upper left 1.267166

% assem4, 3x3, vac 1 side, k =  0.740383




