%> @file  test_erme.m
%> @brief Demo of ERME and solvers.

clear

% Include the MTP source directory
path(path, '../')

% ==============================================================================
% Server (here, on-the-fly)
% ==============================================================================

so = 1; ao = 3;

inputrf = Input();
put(inputrf, 'number_groups',         2);
put(inputrf, 'inner_tolerance',       1e-9);
put(inputrf, 'inner_max_iters',       100);
put(inputrf, 'outer_tolerance',       1e-9);
put(inputrf, 'outer_max_iters',       10);
put(inputrf, 'inner_print_out',       0);
put(inputrf, 'rf_print_out',          0);
put(inputrf, 'quad_order',            2);
put(inputrf, 'rf_max_order_space',       so);
put(inputrf, 'rf_max_order_azimuth',     ao);
put(inputrf, 'rf_max_order_polar',        0);
put(inputrf, 'rf_order_space',       so);
put(inputrf, 'rf_order_azimuth',     ao);
put(inputrf, 'rf_order_polar',        0);
put(inputrf, 'rf_number_nodes', 1);
put(inputrf, 'owner',           'rf');

% Materials
mat = mat_2g_4m();

% Assemblies
[a1,a2,a3,a4] = two_group_benchmark_assemblies();
assemblies = {a4};
server = ResponseServer(inputrf,mat,assemblies);

% ==============================================================================
% ERME setup
% ==============================================================================
input = Input();
put(input, 'inner_tolerance',       1e-10);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-9);
put(input, 'outer_max_iters',       20);
put(input, 'rf_order_space',   so); % Eventually, these will be the
put(input, 'rf_order_azimuth', ao); % orders extracted from a DB
put(input, 'rf_order_polar',    0); 
put(input, 'number_groups',     2);
put(input, 'dimension',         2);
put(input, 'owner',           'erme');

% Define problems.
number_problems = 7;
prob{1}.bc  = {'reflect', 'reflect', 'reflect', 'reflect'};
prob{2}.bc  = {'vacuum',  'reflect', 'reflect', 'reflect'};
prob{3}.bc  = {'reflect', 'vacuum',  'vacuum',  'reflect'};
prob{4}.bc  = {'reflect', 'reflect', 'reflect', 'reflect'};
prob{5}.bc  = {'vacuum',  'reflect', 'reflect', 'reflect'};
prob{6}.bc  = {'reflect', 'vacuum',  'vacuum',  'reflect'};
prob{7}.bc  = {'reflect', 'reflect', 'reflect', 'reflect'};
prob{1}.map = [4 4 4 ; 4 4 4 ; 4 4 4];
prob{2}.map = [4 4 4 ; 4 4 4 ; 4 4 4]/4;
prob{3}.map = [4 4 4 ; 4 4 4 ; 4 4 4];
prob{4}.map = [1 1 1 ; 1 1 1 ; 1 1 1];
prob{5}.map = [1 1 1 ; 1 1 1 ; 1 1 1];
prob{6}.map = [1 1 1 ; 1 1 1 ; 1 1 1];
prob{7}.map = [3 2 1 2 1 ; 2 1 2 1 2 ; 1 2 1 2 1 ; 2 1 2 1 2 ; 1 2 1 2 1];
k = zeros(number_problems, 1);
for i = 2
    disp([' two group benchmark -- erme -- problem ',num2str(i)])
    put(input, 'bc_left',   prob{i}.bc{1});
    put(input, 'bc_right',  prob{i}.bc{2});
    put(input, 'bc_top',    prob{i}.bc{3});
    put(input, 'bc_bottom', prob{i}.bc{4});
    problem = ERME(input, server, prob{i}.map);
    solver = ERME_Picard(input, problem); 
    solve(solver);   
end




