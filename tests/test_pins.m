%> @file  test_pins.m
%> @brief Construct and return simple pin meshes.
function [pin1, pin2, pin3, pin4] = test_pins(verify)

if nargin == 0
    verify = 0;
end

% Note, the material id's corresponsd to those in mat_2g_4m.

% Shared pin properties
Mesh(); % @todo Without this, it complains about Mesh.  why?

pin_pitch = 1.26; radii = 0.54; number = 8;
% Pin 1 - Fuel 1
matid = [2 1];  
pin1 = PinCell(pin_pitch, radii, matid); meshify(pin1, number);
% Pin 2 - Fuel 2
matid = [4 1];  pin2 = PinCell(pin_pitch, radii, matid); meshify(pin2, number);
% Pin 3 - MOD
matid = [  1];  pin3 = PinCell(pin_pitch,    [], matid); meshify(pin3, number);
% Pin 4 - Pure Fuel 2
matid = [  2];  pin4 = PinCell(pin_pitch,    [], matid); meshify(pin4, number);

% Optional verification of pins.
if verify == 1
    % Setup any data needed for evaluation statements
    map = mesh_map(pin2, 'MATERIAL');
    % Define evaluation statements
    tests = {'pitch(pin1)==1.26', 'dx(pin2, 1)==0.1575', 'map(1,1)==1', ...
        'map(4,4)==4','number_cells(pin4)==64','number_cells_x(pin4)==8'};
    run_tests(tests);
end

end

