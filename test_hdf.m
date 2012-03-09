% test of hdf5 writing
%
% simple response function database of the form
%
% /assembly/k/R
% /assembly/k/F  etc.

R1 = rand(10);
F1 = rand(10, 1);
A1 = rand(10, 1);
R2 = rand(10);
F2 = rand(10, 1);
A2 = rand(10, 1);

file = 'sample.h5';

% Number of keffs
h5create(file,'/keffs', [3 1]);

% Orders for all assemblies
h5writeatt(file, '/', 'dimension',   2);
h5writeatt(file, '/', 'solver', 'mtp');
h5writeatt(file, '/', 'max_order_space',   1);
h5writeatt(file, '/', 'max_order_azimuth', 2);
h5writeatt(file, '/', 'max_order_polar',   0);
h5writeatt(file, '/', 'number_keff', 3);
h5write(file, '/keffs', [0.99 1.00 1.01]');

h5create(file,'/a1/k1/R', [10 10]);
h5writeatt(file, '/a1', 'description', 'c5g7 UO2 assembly');
h5writeatt(file, '/a1', 'width_x', 20.0);
h5writeatt(file, '/a1', 'width_y', 20.0);
h5writeatt(file, '/a1', 'width_z',  1.0);
h5write(file, '/a1/k1/R', R1);

h5create(file,'/a1/k1/F', [10 1]);
h5write(file, '/a1/k1/F', F1);

h5create(file,'/a1/k1/A', [10 1]);
h5write(file, '/a1/k1/A', A1);




