function plot_vtk(array, filename)
%  savevtk Save a 2-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.

nx = length(array(:,1));
ny = length(array(1,:));
nz = 1;

fid = fopen(filename, 'wt');
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Comment goes here\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, '\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
%fprintf(fid, '\n');
fprintf(fid, 'ORIGIN     0.000   0.000   0.000\n');
fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
f = reshape(array, nx*ny, 1);

fprintf(fid, 'SCALARS data double\n' );
fprintf(fid, 'LOOKUP_TABLE default\n');

for k=1:nz
    for j=1:ny
        for i=1:nx
            fprintf(fid, '%d ', f(i+(j-1)*nx));
        end
        fprintf(fid, '\n');
    end
end


fclose(fid);
return
end