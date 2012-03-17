%> @file  tests/regression/test_rf_db.m
%> @brief Test of driver and db functions.
function test_rf_db(tester, mat, assemblies)


% Run the driver
inputrf = Input();
put(inputrf, 'number_groups',         2);
put(inputrf, 'inner_tolerance',       1e-12);
put(inputrf, 'inner_max_iters',       100);
put(inputrf, 'outer_tolerance',       1e-12);
put(inputrf, 'outer_max_iters',       10);
put(inputrf, 'inner_print_out',       0);
put(inputrf, 'rf_print_out',          0);
put(inputrf, 'quad_order',            2);
put(inputrf, 'dimension',             2);
put(inputrf, 'rf_max_order_space',    1);
put(inputrf, 'rf_max_order_azimuth',  1);
put(inputrf, 'rf_max_order_polar',    0);
put(inputrf, 'rf_number_nodes',       4);
put(inputrf, 'rf_keff_vector',        [0.9; 1.0]);
put(inputrf, 'rf_db_name',            'test_db.h5');
driver = ResponseDriver(inputrf,mat,assemblies);
driver.run();

%Get all the response
[R, F, A, L] = driver.get_all_responses();
R1(:, :, :) = R(:, :, 2, :);
F1(:, :)    = F(:, 2, :);
A1(:, :)    = A(:, 2, :);
L1(:, :, :) = L(:, :, 2, :);

% Open a DB
inputdb = Input();
put(inputdb, 'rf_db_name',        'test_db.h5');
put(inputdb, 'rf_order_space',    1); % Eventually, these will be the
put(inputdb, 'rf_order_azimuth',  1); % orders extracted from a DB
put(inputdb, 'rf_order_polar',    0); 
put(inputdb, 'rf_interp_method',  'spline'); 
db = ResponseDB(inputdb);

% Read the responses
read_response(db);

% Get responses at known values of keff
[R2, F2, A2, L2] = get_responses(db, 1.0);

% Define evaluation statements and test them.
tests = {'tester.almost_equal(R1(1,1,1), 0.150900409597154)', ...
         'tester.almost_equal(F1(1,2),   0.033936809520236)', ...
         'tester.almost_equal(A1(1,3),   0.035917412550710)', ...
         'tester.almost_equal(L1(1,1,4), 0.229660401229798)', ...
         'tester.almost_equal(norm(reshape(R1, prod(size(R1)), 1)-reshape(R2, prod(size(R2)), 1)), 0.0)', ...
         'tester.almost_equal(norm(F1-F2), 0.0)', ...
         'tester.almost_equal(norm(A1-A2), 0.0)', ...
         'tester.almost_equal(norm(reshape(L1, prod(size(L1)), 1)-reshape(L2, prod(size(L2)), 1)), 0.0)'};
tester.run_tests(tests);





