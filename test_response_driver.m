% test of response driver

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'inner_tolerance',       1e-9);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-9);
put(input, 'outer_max_iters',       10);
put(input, 'quad_order',            2);
put(input, 'rf_max_order_space',        2);
put(input, 'rf_max_order_azimuth',      2);
put(input, 'rf_max_order_polar',        0);
put(input, 'rf_k_vector', 0.9009968666);

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 10;
% Pin 1 - Fuel 1
matid = [2 1];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2 - Fuel 2
matid = [4 1];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3 - MOD
matid = [  2];  pin3 = PinCell(pitch,    [], matid); meshify(pin3, number);

% ==============================================================================
% ASSEMBLIES
% ==============================================================================

% Assembly 1 
pin_map1 = [1 1 1 
            1 2 1
            1 1 1];     
mesh1 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh1);

% Assembly 2 
pin_map1 = [1 1 1 
            1 3 1
            1 1 1];     
mesh2 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh2);

% ==============================================================================
% RESPONSE FUNCTION DRIVER
% ==============================================================================
driver = ResponseDriver(input, mat, {mesh1});
run(driver);
