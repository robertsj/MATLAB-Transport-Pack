function mesh = test_mesh(n)

if nargin == 0
    n = 1;
end

% ------------------------------------------------------------------------------
% Simple 1 group, 2 material
% ------------------------------------------------------------------------------ 
if n == 1
    a=0.5;
    % Simple uniform box.
    xcm = [ 0.0  10.0 ]*a;
    xfm = [   10   ]*a;
    ycm = [ 0.0  10.0 ]*a;
    yfm = [   10   ]*a;    
    % Material map
    mat_map = [1];
    % Make the mesh.
    mesh = Mesh2D(xfm, yfm, xcm, ycm, mat_map);
    
elseif n == 2
    
    xcm    = [ 0.0  10.0  25.0];
    xfm    = [    10    20    ];
    ycm    = [ 0.0  15.0  30.0];
    yfm    = [    20    40    ]/2;

    mat_map = [ 1  3     % ^
                2  2 ];  % | y   x -->

    mesh = Mesh2D(xfm, yfm, xcm, ycm, mat_map);

    DBC.Assert('dx(mesh,  1) == 1.0')
    DBC.Assert('dy(mesh, 21) == 15/40')
    DBC.Assert('number_cells_x(mesh) == 30')
    DBC.Assert('number_cells_y(mesh) == 60')

    DBC.Assert('mat_map_fine( 1,  1) == 2')
    DBC.Assert('mat_map_fine(11,  1) == 2')
    DBC.Assert('mat_map_fine( 1, 21) == 1')
    DBC.Assert('mat_map_fine(11, 21) == 3')
    
elseif n == 3
    
    xcm    = [ 0.0  20     40]/4;
    xfm    = [    5 5  ];
    ycm    = [ 0.0  20    40]/4;
    yfm    = [    5 5  ];

    mat_map = [ 2 2
                2 2];  % | y   x -->

    mesh = Mesh2D(xfm, yfm, xcm, ycm, mat_map);

end