%> @file  test_pins.m
%> @brief Construct and return simple pin meshes.
function [pin1, pin2, pin3, pin4] = two_group_benchmark_pins(verify)

if nargin == 0
    verify = 0;
end

% Note, the material id's corresponsd to those in mat_2g_4m.


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

