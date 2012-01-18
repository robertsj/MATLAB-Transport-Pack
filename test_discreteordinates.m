function test_discreteordinates
% This is a test case for DiscreteOrdinates that demonstrates basic use of the
% class and solver options.  Note, a lot of the basic usage is similar to that
% of Denovo in its Python form.

% Create the mesh.
mesh = test_mesh();

% Create the materials.
mat = test_materials();

% Create an external source.
source = test_source();

% Create the manager.
manager = DiscreteOrdinates();

end % end test_discreteordinates
