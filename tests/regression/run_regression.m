%> @file  tests/regression/run_regression.m
%> @brief Run a few test cases to ensure all is well.

clear all

% Include base test path.
path(path, '../')

% Initialize the test driver.
tester = TestDriver();

% Test materials (2 group, 4 material)
mat = test_materials(tester);

% Test pins 
pins = test_pins(tester);

% Test assemblies
assemblies = test_assemblies(tester, pins);

% Test response driver and db
test_rf_db(tester, mat, assemblies);

% Total
tester.results();