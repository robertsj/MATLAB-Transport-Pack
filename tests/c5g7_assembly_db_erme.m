%> @file  test_erme.m
%> @brief Demo of ERME and solvers.

clear

% Include the MTP source directory
path(path, '../')

tic

input = Input();

so = [0 1 2 4];
ao = [0 2 4];
po = [0 2];

for s = 4
    for a = 3
        for p = 2

            
            disp([' so = ',num2str(s),' ao = ',num2str(a),' po = ',num2str(p)])
% ==============================================================================
% Server
% ==============================================================================
put(input, 'rf_db_name', 'c5g7_assemblies.h5');
put(input, 'rf_interp_method', 'spline');
put(input, 'rf_order_space',        so(s)); % Eventually, these will be the
put(input, 'rf_order_azimuth',      ao(a)); % orders extracted from a DB
put(input, 'rf_order_polar',        po(p)); 
server = ResponseDB(input);
read_response(server);

% ==============================================================================
% ERME setup
% ==============================================================================
put(input, 'inner_tolerance',       1e-10);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-9);
put(input, 'outer_max_iters',       20);
put(input, 'number_groups',         7);
put(input, 'dimension',             2);
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'vacuum');

map = [ 1 2 3
        2 1 3
        3 3 3 ];
problem = ERME(input, server, map);
solver  = ERME_Picard(input, problem); 
solve(solver);   

toc

        end
    end
end



