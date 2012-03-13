% test of response driver
clear all
% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         2);
put(input, 'inner_tolerance',       1e-10);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-10);
put(input, 'outer_max_iters',       10);
put(input, 'quad_order',            2);
put(input, 'rf_max_order_space',       13);
put(input, 'rf_max_order_azimuth',      3);
put(input, 'rf_max_order_polar',        0);
kv = linspace(0.5, 1.4, 19);
put(input, 'rf_k_vector',   kv');
put(input, 'rf_number_nodes', 4);
put(input, 'rf_db_name', 'two_group_3x3.h5');

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(2);

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 8;
% Pin 1 - Fuel 1
matid = [2 1];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2 - Fuel 2
matid = [4 1];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3 - MOD
matid = [  1];  pin3 = PinCell(pitch,    [], matid); meshify(pin3, number);
% Pin 4 - Pure Fuel 2
matid = [  2];  pin4 = PinCell(pitch,    [], matid); meshify(pin4, number);

% ==============================================================================
% ASSEMBLIES
% ==============================================================================

% Assembly 1 (kinf =  1.319782)  1.319782 0.891295 1.297046 | 1.278287
pin_map1 = [1 1 1 
            1 1 1
            1 1 1];     
mesh1 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh1);

% Assembly 2 (kinf =  0.891295)
pin_map1 = [1 1 1 
            1 2 1
            1 1 1];     
mesh2 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh2);

% Assembly 3 (kinf =  1.297046)
pin_map1 = [1 1 1 
            1 3 1
            1 1 1];     
mesh3 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh3);

% Assembly 4 (kinf =   1.278287)
pin_map1 = [4 4 4
            4 4 4
            4 4 4];     
mesh4 = Assembly({pin1, pin2, pin3, pin4}, pin_map1); 
meshify(mesh4);

% ==============================================================================
% RESPONSE FUNCTION DRIVER
% ==============================================================================
mesh_array = {mesh1, mesh2, mesh3, mesh4};
driver = ResponseDriver(input, mat, mesh_array);
run(driver);
[R, F, A, L] = get_responses(driver);
% RR=R(:,:,1,1);
% elements = [1];
% number_elements = 1;        
% M = Connect(input, elements, number_elements);
% [v, e]=eigs(M*RR); 
% J = v(:, 1);
% J = sign(J(1))*J;
% FF = F(:, 1, 1);
% AA = A(:, 1, 1); 
% e = eigs(M*RR, 4, 'LR')
% LL = L(:, :, 1, 1);
% leak = LL'*J;
%  (FF'*J )/(AA'*J + leak(2) + leak(4))
%  
rf_db = ResponseDB(input);
initialize_write(rf_db);

for i = 1:length(mesh_array)
   write_response(rf_db, i, ['assembly',num2str(i)], ...
       R(:,:,:,i), F(:,:,i), A(:,:,i), L(:,:,:,i));
end

read_response(rf_db);

