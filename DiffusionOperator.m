%> @file  DiffusionOperator.m
%> @brief DiffusionOperator class definition.
% ==============================================================================
%> @brief Mesh-centered finite difference diffusion operator.
% ==============================================================================
classdef DiffusionOperator < handle
    
    properties
        d_input
        d_mat
        d_mesh
        d_mg_operator
        d_1g_operators
        d_built = 0;
        d_nx
        d_ny
        d_nz
        d_n
        d_albedo
    end
    
    methods
        
        % ======================================================================
        %> @brief  Class constructor
        %> @param  input    Input database.
        %> @param  mat      Material database.
        %> @param  mesh     Cartesian mesh.
        %> @return          Instance of the DiffusionOperator class.
        % ======================================================================
        function this = DiffusionOperator(input, mat, mesh) 
            this.d_input = input;
            this.d_mat   = mat;
            this.d_mesh  = mesh;
            
            this.d_nx = number_cells_x(mesh);
            this.d_ny = number_cells_y(mesh);
            this.d_nz = number_cells_z(mesh);
            this.d_n  = number_cells(mesh);
            
            this.d_albedo = zeros(6, 1);
            if strcmp(input.get('bc_left'), 'reflect')
                this.d_albedo(1) = 1;
            end
            if strcmp(input.get('bc_right'), 'reflect')
                this.d_albedo(2) = 1;
            end
            if strcmp(input.get('bc_bottom'), 'reflect') || mesh.DIM < 2
                this.d_albedo(3) = 1;
            end
            if strcmp(input.get('bc_top'), 'reflect') || mesh.DIM < 2
                this.d_albedo(4) = 1;
            end
            if strcmp(input.get('bc_south'), 'reflect') || mesh.DIM < 3
                this.d_albedo(5) = 1;
            end            
            if strcmp(input.get('bc_north'), 'reflect') || mesh.DIM < 3
                this.d_albedo(6) = 1;
            end            
            
        end 
        
        function M = get_1g_operator(this, g)
            if ~this.d_built
                build_1g_operators(this);
            end
            M = this.d_1g_operators{g};
        end
        
        function M = get_mg_operator(this)
            if ~this.d_built
                build_mg_operator(this);
            end
            M = this.d_mg_operator;  
        end
        
    end
    
    methods (Access = private)
        
        function this = build_mg_operator(this)
            
        end
        
        % S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to 
        % generate an m-by-n sparse matrix such that S(i(k),j(k)) = s(k), 
        % with space allocated for nzmax nonzeros. Vectors i, j, and s are 
        % all the same length. Any elements of s that are zero are ignored, 
        % along with the corresponding values of i and j. Any elements of s 
        % that have duplicate values of i and j are added together.
        
        function this = build_1g_operators(this)
            
            % Initialize cell array of 1g operators.
            this.d_1g_operators = cell(number_groups(this.d_mat), 1);
            
            % Get material mesh map;
            map = reshape(mesh_map(this.d_mesh, 'MATERIAL'), ...
                          number_cells(this.d_mesh), 1);
            
            % Simplify notation.
            mat  = this.d_mat;        
            mesh = this.d_mesh; 
            
            count = 0;
            
            nxyz = [1 this.d_nx
                    1 this.d_ny
                    1 this.d_nz];
            
            for g = 1:number_groups(mat)
                
                % Presize the operator arrays.  Do it better later.
                value   = zeros(this.d_n*(1+1*mesh.DIM), 1);
                index_i = value+1;
                index_j = value+1;
                
                for row = 1:this.d_n
                    
                    % This cell's data
                    cell_dc   = 1 / (3*mat.sigma_t(map(row), g));
                    cell_sigr = mat.sigma_t(map(row), g) - ...
                                mat.sigma_s(map(row), g, g);                   
                    [i, j, k] = matrix_to_indices(this, row);          
                    cell_hxyz = mesh.dx(i)*mesh.dy(j)*mesh.dz(k);
                    
                    bound = [i i j j k k];
                    
                    jo = zeros(6, 1);
                    
                    % leak=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z
                    for leak = 1:6
                        
                        % define (x,y,z) and (-,+) indices
                        xyz_idx   = ceil(leak/2);        % x, y, or z
                        dir_idx   = 2 - mod(leak, 2);    % -=1, +=2 (left/right)
                        
                        neig_idx  = [i j k];             % begin with i,j,k
                        shift_idx = -2*mod(leak, 2) + 1; % shift neig by -1 or +1
                        neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx);
                        
                        % Get neigbor cell data.
                        neig_row = index(this.d_mesh, ...
                            neig_idx(1), neig_idx(2), neig_idx(3));
                        [ii, jj, kk] = matrix_to_indices(this, neig_row);

                        
                        % Compute coupling coefficient
                        if (bound(leak) == nxyz(xyz_idx,dir_idx))
                            
                            dtilde = (2*cell_dc*(1-this.d_albedo(leak))) / ...
                                     (4*cell_dc*(1+this.d_albedo(leak)) + ...
                                      (1-this.d_albedo(leak))*cell_hxyz);
                            
                        else  % not a boundary
                            
                            % Neighbor volume and diffusion coefficient.
                            try
                            neig_hxyz = mesh.dx(ii)*mesh.dy(jj)*mesh.dz(kk);
                            catch
                               a=1; 
                            end
                            neig_dc   = 1 / (3*mat.sigma_t(map(neig_row), g));
                            
                            % Compute dtilde.
                            dtilde = (2*cell_dc*neig_dc) / ...
                                     (neig_hxyz*cell_dc + cell_hxyz*neig_dc);
                            
                            % Add to matrix arrays.
                            count = count + 1;
                            value(count)   = -dtilde / cell_hxyz;
                            index_i(count) = row;
                            index_j(count) = neig_row;
                            
                        end
                        
                        
                        % compute leakage coefficient for target
                        jo(leak) = shift_idx*dtilde;
                        
                    end % leak loop
                    
                    jnet = (jo(2) - jo(1))/mesh.dx(i) + ...
                           (jo(4) - jo(3))/mesh.dy(j) + ...
                           (jo(6) - jo(5))/mesh.dz(k);

                    % Add to matrix arrays.
                    count = count + 1;
                    value(count)   = jnet + cell_sigr;
                    index_i(count) = row;
                    index_j(count) = row;
                    
                    
                end % row loop
                
                this.d_1g_operators{g} = sparse(index_i, index_j, value);

            end % group loop
            
        end
        
        function [i, j, k] = matrix_to_indices(this, row)
            i = mod(row-1, this.d_nx) + 1;
            j = floor(mod(row-1, this.d_nx*this.d_ny)/(this.d_nx)) + 1;
            k = floor(mod(row-1, this.d_nx*this.d_ny*this.d_nz)/(this.d_nx*this.d_ny)) + 1;
        end
           
        
    end
    
end