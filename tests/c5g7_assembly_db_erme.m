%> @file  test_erme.m
%> @brief Demo of ERME and solvers.

clear

% Include the MTP source directory
path(path, '../')

tic

input = Input();

so = [0 1 2 3 4];
ao = [0 1 2 3 4];
po = [0 1 2];

% for s = 1:length(so)
%     for a = 1:length(ao)
%         for p = 1:length(po)
s = 2;
a = 2;
p = 3;
            
            %disp([' so = ',num2str(s),' ao = ',num2str(a),' po = ',num2str(p)])
% ==============================================================================
% Server
% ==============================================================================
put(input, 'rf_db_name',            'c5g7_assemblies.h5');
put(input, 'rf_db_name_mat',        'c5g7_assemblies.mat');
put(input, 'rf_interp_method',      'linear');
put(input, 'rf_order_space',        so(s)); % Eventually, these will be the
put(input, 'rf_order_azimuth',      ao(a)); % orders extracted from a DB
put(input, 'rf_order_polar',        po(p)); 
server = ResponseDB(input);
read_response(server);
%read_response_mat(server);

% ==============================================================================
% ERME setup
% ==============================================================================
put(input, 'inner_tolerance',       1e-10);
put(input, 'inner_max_iters',       2000);
put(input, 'outer_tolerance',       1e-9);
put(input, 'outer_max_iters',       20);
put(input, 'number_groups',         7);
put(input, 'dimension',             2);
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'vacuum');
put(input, 'bc_top',                'reflect');
put(input, 'bc_bottom',             'vacuum');
put(input, 'erme_steffensen',       0);
put(input, 'erme_picard_inner',     'eigs');
put(input, 'slepc_options',         {'-malloc','-malloc_debug','-malloc_dump'});
map = [ 1 2 3
        2 1 3
        3 3 3 ];
problem = ERME(input, server, map);
solver  = ERME_Picard(input, problem); 

t = toc;
format long
solve(solver);   
toc
problem.get_keff()
problem.get_lambda()
% keff2(s, a, p) = problem.get_keff();
% lambda2(s,a,p) = problem.get_lambda();
% time2(s, a, p) = toc - t;
% its2(s, a, p)  = solver.d_number_iterations;

%         end
%     end
% end

kref = 1.187996850234024;


