function q_e = test_source(~)

% Get the test mesh.
mesh = test_mesh();

% Make the source.
q_e = Source(mesh, 3);

% Make two unique, three-group sources.
spectra = [ 1.1  0.3 0.1
            2.2  0.5 0.1
            3.3  0.7 0.1];

placement = [ 1 2
              2 3 ];
          
q_e = set_sources(q_e, spectra, placement);          

if nargin == 1
    figure(1)
    plot_mesh_map(mesh, 'EXTERNALSOURCE');
    figure(2)
    plot_mesh_map(mesh, 'MATERIAL');
end

DBC.Assert('source(q_e, 1, 1, 1) == 0.3/4.0/pi');

end