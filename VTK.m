%> @file  VTK.m
%> @brief Utilities to print data in VTK format for viewing in Visit etc.
classdef VTK
    
    methods (Static)
        
        function savevtk(state, ng, mesh, filename)
            %  savevtk Save a 2-D scalar array in VTK format.
            %  savevtk(array, filename) saves a 3-D array of any size to
            %  filename in VTK format.

            nx = number_cells_x(mesh);
            ny = number_cells_y(mesh);
            nz = number_cells_z(mesh);
            
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
            for g = 1:ng
                f = flux(state, g);
                %fprintf(fid, '\n');
                
                fprintf(fid, 'SCALARS group%1i double\n', g);
                fprintf(fid, 'LOOKUP_TABLE default\n');
                %fprintf(fid, '\n');
                for k=1:nz
                    for j=1:ny
                        for i=1:nx
                            fprintf(fid, '%d ', f(index(mesh, i, j, k)));
                        end
                        fprintf(fid, '\n');
                    end
                end
            end
            
            fclose(fid);
            return
        end
        
        function savevtk3d(array, filename)
            %  savevtk Save a 3-D scalar array in VTK format.
            %  savevtk(array, filename) saves a 3-D array of any size to
            %  filename in VTK format.
            [nx, ny, nz] = size(array);
            fid = fopen(filename, 'wt');
            fprintf(fid, '# vtk DataFile Version 2.0\n');
            fprintf(fid, 'Comment goes here\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, '\n');
            fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
            fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
            fprintf(fid, '\n');
            fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
            fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
            fprintf(fid, '\n');
            fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
            fprintf(fid, 'SCALARS scalars double\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '\n');
            for a=1:nz
                for b=1:ny
                    for c=1:nx
                        fprintf(fid, '%d ', array(c,b,a));
                    end
                    fprintf(fid, '\n');
                end
            end
            fclose(fid);
            return
        end
            
        function savevtkvector(X, Y, Z, filename)
            %  savevtkvector Save a 3-D vector array in VTK format
            %  savevtkvector(X,Y,Z,filename) saves a 3-D vector of any size to
            %  filename in VTK format. X, Y and Z should be arrays of the same
            %  size, each storing speeds in the a single Cartesian directions.
            if ((size(X) ~= size(Y)) | (size(X) ~= size(Z)))
                fprint('Error: velocity arrays of unequal size\n'); return;
            end
            [nx, ny, nz] = size(X);
            fid = fopen(filename, 'wt');
            fprintf(fid, '# vtk DataFile Version 2.0\n');
            fprintf(fid, 'Comment goes here\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, '\n');
            fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
            fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
            fprintf(fid, '\n');
            fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
            fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
            fprintf(fid, '\n');
            fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
            fprintf(fid, 'VECTORS vectors double\n');
            fprintf(fid, '\n');
            for a=1:nz
                for b=1:ny
                    for c=1:nx
                        fprintf(fid, '%f ', X(c,b,a));
                        fprintf(fid, '%f ', Y(c,b,a));
                        fprintf(fid, '%f ', Z(c,b,a));
                    end
                    fprintf(fid, '\n');
                end
            end
            fclose(fid);
            return
        end
        
    end
    
end