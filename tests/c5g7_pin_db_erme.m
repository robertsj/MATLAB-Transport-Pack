%> @file  c5g7_pin_db_erme.m
%> @brief Generate response database for C5G7 pin cells.
%

clear

% Include the MTP source directory
path(path, '../')

% ==============================================================================
% Create a response DB
% ==============================================================================

inputrf = Input();
put(inputrf, 'number_groups',         7);
put(inputrf, 'inner_tolerance',       1e-12);
put(inputrf, 'inner_max_iters',       100);
put(inputrf, 'outer_tolerance',       1e-9);
put(inputrf, 'outer_max_iters',       10);
put(inputrf, 'inner_print_out',       0);
put(inputrf, 'rf_print_out',          1);
put(inputrf, 'quad_order',            18);
put(inputrf, 'rf_max_order_space',    4);
put(inputrf, 'rf_max_order_azimuth',  4);
put(inputrf, 'rf_max_order_polar',    1);
put(inputrf, 'rf_keff_vector',           [0.9 1.0 1.10 1.20 1.3]');
%put(inputrf, 'rf_k_vector',           [1.0]');

put(inputrf, 'rf_number_nodes',       7);
put(inputrf, 'rf_db_name', 'c5g7_pins.h5');

% ==============================================================================
% MATERIALS
% ==============================================================================
mat = C5G7_materials(0);

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 13; %52;
% Pin 1 - UO2 
matid = [1 7];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2 - 4.3% MOX
matid = [2 7];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3 - 7.0% MOX
matid = [3 7];  pin3 = PinCell(pitch, radii, matid); meshify(pin3, number);
% Pin 4 - 8.7% MOX
matid = [4 7];  pin4 = PinCell(pitch, radii, matid); meshify(pin4, number);
% Pin 5 - Guide Tube
matid = [5 7];  pin5 = PinCell(pitch, radii, matid); meshify(pin5, number);
% Pin 6 - Fission Chamber
matid = [6 7];  pin6 = PinCell(pitch, radii, matid); meshify(pin6, number);
% Pin 7 - Moderator
matid = [7  ];  pin7 = PinCell(pitch, [   ], matid); meshify(pin7, number);

%plot_mesh_map(pin1, 'MATERIAL')

%return

% ==============================================================================
% RESPONSE FUNCTION DRIVER
% ==============================================================================
mesh_array = {pin1, pin2, pin3, pin4, pin5, pin6, pin7};
driver = ResponseDriver(inputrf, mat, mesh_array);
run(driver);
[R, F, A, L] = get_responses(driver);


% ==============================================================================
% RESPONSE FUNCTION DATABASE
% ==============================================================================

rf_db = ResponseDB(inputrf);
initialize_write(rf_db);
for i = 1:length(mesh_array)
   write_response(rf_db, i, ['pin',num2str(i)], ...
       R(:,:,:,i), F(:,:,i), A(:,:,i), L(:,:,:,i));
end





