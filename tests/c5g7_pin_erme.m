%> @file  c5g7_pin_erme.m
%> @brief Solve c5g7 via ERME and pin cell responses

clear
% Include the MTP source directory
path(path, '../')

% ==============================================================================
% Open DB
% ==============================================================================
input = Input();
put(input, 'rf_db_name',         'c5g7_pins.h5');
put(input, 'rf_order_space',    1); % Eventually, these will be the
put(input, 'rf_order_azimuth',  2); % orders extracted from a DB
put(input, 'rf_order_polar',    1); 
put(input, 'rf_interp_method',  'spline'); 
server = ResponseDB(input);
read_response(server);

% ==============================================================================
% Setup ERME
% ==============================================================================
put(input, 'bc_left',           'reflect');
put(input, 'bc_right',          'reflect');
put(input, 'bc_bottom',         'reflect');
put(input, 'bc_top',            'reflect');
put(input, 'inner_tolerance',   1e-10);
put(input, 'inner_max_iters',   100);
put(input, 'outer_tolerance',   1e-9);
put(input, 'outer_max_iters',   20);
put(input, 'number_groups',     7);
put(input, 'dimension',         2);

% ==============================================================================
% ASSEMBLIES
% ==============================================================================

G = 5; F = 6;

% Assembly 1 - UO2
pin_map1 = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 G 1 1 G 1 1 G 1 1 1 1 1 
            1 1 1 G 1 1 1 1 1 1 1 1 1 G 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 G 1 1 G 1 1 G 1 1 G 1 1 G 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 G 1 1 G 1 1 F 1 1 G 1 1 G 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 G 1 1 G 1 1 G 1 1 G 1 1 G 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 G 1 1 1 1 1 1 1 1 1 G 1 1 1 
            1 1 1 1 1 G 1 1 G 1 1 G 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];     
% Assembly 2 - MOX
pin_map2 = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
            2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2
            2 3 3 3 3 G 3 3 G 3 3 G 3 3 3 3 2
            2 3 3 G 3 4 4 4 4 4 4 4 3 G 3 3 2
            2 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 2
            2 3 G 4 4 G 4 4 G 4 4 G 4 4 G 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 G 4 4 G 4 4 F 4 4 G 4 4 G 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 3 4 4 4 4 4 4 4 4 4 4 4 3 3 2
            2 3 G 4 4 G 4 4 G 4 4 G 4 4 G 3 2
            2 3 3 3 4 4 4 4 4 4 4 4 4 3 3 3 2
            2 3 3 G 3 4 4 4 4 4 4 4 3 G 3 3 2
            2 3 3 3 3 G 3 3 G 3 3 G 3 3 3 3 2
            2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2
            2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
% Assembly 3 - Moderator
pin_map3 = [7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
            7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7];
        
elements = [1 1 1 1 1 1 1 1 
            1 1 1 1 1 1 1 1  
            1 1 1 1 1 G 1 1 
            1 1 1 G 1 1 1 1 
            1 1 1 1 1 1 1 1 
            1 1 G 1 1 G 1 1  
            1 1 1 1 1 1 1 1  
            1 1 1 1 1 1 1 1];
        
problem = ERME(input, server, elements);
solver = ERME_Picard(input, problem);
tic, solve(solver); toc


