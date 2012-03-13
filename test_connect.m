%> @file  test_connect.m
%> @brief Demo of 1/2/3-d connectivity matrix.

% ==============================================================================
% Input
% =============================================================================
clear 

input = Input();
put(input, 'bc_left',           'reflect');
put(input, 'bc_right',          'reflect');
put(input, 'bc_top',            'reflect');
put(input, 'bc_bottom',         'reflect');
put(input, 'bc_north',          'reflect');
put(input, 'bc_south',          'reflect');
put(input, 'rf_order_space',    3); % Eventually, these will be the
put(input, 'rf_order_azimuth',  1); 
put(input, 'rf_order_polar',    1); % orders extracted from a DB
put(input, 'number_groups',     2);

% 1-d
put(input, 'dimension',     1);
e1 = [1 2 3];
b1 = Connect(input, e1);
M1 = build(b1);
figure(1)
plot_connect(b1)

% 2-d
put(input, 'dimension',     2);
e2 = [1 1
      2 1]; 
e2 = flipud(e2)';
b2 = Connect(input, e2);
M2 = build(b2);
figure(2)
plot_connect(b2)


% 3-d 
put(input, 'dimension',     3);
e3(:, :, 1) = [1 1
               1 1]; 
e3(:, :, 2) = [1 1
               1 1]; 
b3 = Connect(input, e3);
M3 = build(b3);
figure(3)
plot_connect(b3)