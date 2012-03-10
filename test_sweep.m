



% Test input.
input = Input();

% Get mesh and materials
mesh = test_mesh(1);
mat  = test_materials(2);

% Quadrature
quadrature = LevelSymmetric(16);

% Boundary
boundary = BoundaryMesh(input, mesh, quadrature);
initialize(boundary, 1);

% Make DD equation
equation = DD2D(mesh, mat, quadrature);
setup_group(equation, 1);

% Make sweeper
sweeper = Sweep2D_mod(input, mesh, mat, quadrature, boundary, equation);
sweeper.d_kernel = 0;
% Make a sweep source
source  = ones(number_cells(mesh), 1);

n = 2;
t  = zeros(n, 1);
for i = 1:n
    tic
    phi     = sweep(sweeper, source, 1);
    t(i) = toc;
end
tmean = mean(t);
tstd  = std(t);
disp([' tmean = ',num2str(tmean), ' +/- ', num2str(tstd)])
sweeper.d_kernel = 1;
t  = zeros(n, 1);
for i = 1:n
    tic
    phi     = sweep(sweeper, source, 1);
    t(i) = toc;
end
tmean = mean(t);
tstd  = std(t);
disp([' tmean = ',num2str(tmean), ' +/- ', num2str(tstd)])


%plot_flux(mesh, phi)
