%> @file  test_erme.m
%> @brief Demo of ERME and solvers.

% ==============================================================================
% Server (here, a DB)
% ==============================================================================
clear
db_input = Input();
put(db_input, 'rf_db_name', 'two_group_3x3.h5');
%put(db_input, 'rf_db_name', 'homog.h5');
rf_db = ResponseDB(db_input);
read_response(rf_db);

% ==============================================================================
% ERME setup
% ==============================================================================
input = Input();
put(input, 'bc_left',           'vacuum');
put(input, 'bc_right',          'reflect');
put(input, 'bc_bottom',         'reflect');
put(input, 'bc_top',            'reflect');
put(input, 'inner_tolerance',       1e-10);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-10);
put(input, 'outer_max_iters',       10);
put(input, 'rf_order_space',   13); % Eventually, these will be the
put(input, 'rf_order_polar',    3); % orders extracted from a DB
put(input, 'rf_order_azimuth',  0); 
put(input, 'number_groups',     2);
put(input, 'dimension',         2);
elements = [4 4 4; 4 4 4; 4 4 4];     
problem = ERME(input, rf_db, elements);
solver = ERME_Picard(input, problem);
tic, solve(solver); toc
% 1.319782 0.891295 1.297046 | 1.278287 <-- 1x1 reflective
% 2g, assem1, vac 1 side, k = 0.551020, 3x3 assem1 k =  0.781512
% 3x3 assem1 + assem2 upper left 1.267166

% assem4, 3x3, vac 1 side, k =  0.740383




