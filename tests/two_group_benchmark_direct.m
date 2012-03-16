%> @file  two_group_benchmark_direct.m
%> @brief Solve several small assembly problems that have known solutions.

% Include the MTP source directory
path(path, '../')

clear all

% Common input
input = Input();
put(input, 'number_groups',         2);
put(input, 'eigen_tolerance',       1e-9);
put(input, 'eigen_max_iters',       2000);
put(input, 'inner_tolerance',       1e-9);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-9);
put(input, 'outer_max_iters',       10);
put(input, 'inner_solver',          'SI');
put(input, 'print_out',             0);
% Quadrature (2nd order Quadruple Range = 2 angles/quadrant)
quadrature  = QuadrupleRange(2);
% Materials
mat = mat_2g_4m();
% Assemblies
[a1,a2,a3,a4] = two_group_benchmark_assemblies();
assemblies = {a1,a2,a3,a4};
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
prob{2}.map = [4 4 4 ; 4 4 4 ; 4 4 4];
prob{3}.map = [4 4 4 ; 4 4 4 ; 4 4 4];
prob{4}.map = [1 1 1 ; 1 1 1 ; 1 1 1];
prob{5}.map = [1 1 1 ; 1 1 1 ; 1 1 1];
prob{6}.map = [1 1 1 ; 1 1 1 ; 1 1 1];
prob{7}.map = [3 2 1 2 1 ; 2 1 2 1 2 ; 1 2 1 2 1 ; 2 1 2 1 2 ; 1 2 1 2 1];
k = zeros(number_problems, 1);
for i = 1:number_problems
    core = Core(assemblies, prob{i}.map);
    meshify(core);
    put(input, 'bc_left',   prob{i}.bc{1});
    put(input, 'bc_right',  prob{i}.bc{2});
    put(input, 'bc_top',    prob{i}.bc{3});
    put(input, 'bc_bottom', prob{i}.bc{4});
    state    = State(input, core);
    boundary = BoundaryMesh(input, core, quadrature);
    q_e      = Source(core, 2);
    q_f      = FissionSource(state, core, mat);
    initialize(q_f);
    solver= Eigensolver(input,      ...
                        state,    	...
                        boundary,   ...
                        core,     	...
                        mat,        ...
                        quadrature, ...
                        q_e,        ...
                        q_f);
    solve(solver);    
    prob{i}.state = state;
    prob{i}.core  = core;
    fprintf(1,' problem = %d, keff = %20.14f \n', i, eigenvalue(state));
end

save('two_benchmark_direct.mat', 'prob');

%     1.27828651169820 
%     0.74038277197922 
%     0.50656311677175 
%     1.31978206812546 
%     0.78151161472123 
%     0.53121197138664 
%     1.08741028492215 

