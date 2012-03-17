%> @file  ResponseServer.m
%> @brief Produces responses on-the-fly.
% ==============================================================================
%> @brief ResponseServer class definition.
%
%> 
% ==============================================================================
classdef ResponseServer < ResponseServerBase

    properties (Access = private)

        %> Material database
        d_mat
        %> Cell array of mesh objects
        d_mesh_array
        %> Cell array of local solver objects
        d_solver
        %> Cell array of local solver fission sources
        d_fission_source
        %> Angular mesh object
        d_quadrature
        %> Boundary object for all problems (*assumes* identical meshing)
        d_boundary
        %> Incident condition (*assumes* just one incident side)
        d_bc
        %> State object for all problems
        d_state
        
        
        d_keff_current = -1;
        d_keff_last = -1;

        %> Working array for fission rates
        d_fiss
        %> Working array for absorption rates
        d_abso
        %> Working array for leakage rates
        d_leak
        %> Precomputed fine mesh fission cross section vector
        d_fission
        %> Precomputed fine mesh absorption cross section vector
        d_absorption
        %> Precomputed fine mesh volume vector
        d_volume
        %> Working cell array of arrays of boundary responses
        d_coef

        %> Spatial basis
        d_basis_space
        %> Polar angle basis        
        d_basis_polar
        %> Azimuthal angle basis        
        d_basis_azimuth
        %> Number of spatial cells   
        d_number_space
        %> Number of polar angles   
        d_number_polar
        %> Number of azimuthal angles
        d_number_azimuth        

    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %> @param input         User input database
        %> @param mat           Material database
        %> @param mesh_array    Cell array of mesh objects
        % ======================================================================
        function this = ResponseServer(input, mat, mesh_array)
            
            % Call base class
            this = this@ResponseServerBase(input);
            
            % Verify input and materials have same group count.
            DBC.Require('get(input, ''number_groups'') == number_groups(mat)');
            % Assert we got at least one mesh.
            DBC.Require('~isempty(mesh_array)');
            % We need the quadrature order.
            DBC.Require('get(input, ''quad_order'')');
            
            % We're limiting studies to 4-way symmetric nodes.
            put(input, 'bc_left',               'response');
            put(input, 'bc_right',              'vacuum');
            put(input, 'bc_top',                'vacuum');
            put(input, 'bc_bottom',             'vacuum');
            
            this.d_input = input;
            this.d_mat   = mat;
            this.d_mesh_array = mesh_array;
            
            % Setup the quadrature.  Could defaulted to 2, but be explicit!
            qo = get(input, 'quad_order');
            this.d_quadrature  = QuadrupleRange(qo);
            
            % Get the orders we'll use.
            this.d_number_groups = get(this.d_input, 'number_groups');
            this.d_order_space   = get(this.d_input, 'rf_order_space');
            this.d_order_azimuth = get(this.d_input, 'rf_order_azimuth');
            this.d_order_polar   = get(this.d_input, 'rf_order_polar');            
            
            % I'm not sure how to generalize for MOC at the moment.  For meshes,
            % it's assumed all nodes have same spacing.  This could be a future
            % extension.
            this.d_boundary = ...
                BoundaryMesh(input, mesh_array{1}, this.d_quadrature);
            
            % Get the left boundary condition, the only one used for now.
            this.d_bc = get_bc(this.d_boundary, Mesh.LEFT);
            
            % Create the state vector based on one of possibly several
            % identically-spaced meshes.
            this.d_number_nodes = length(mesh_array);
            this.d_state  = State(input, mesh_array{1});
            
            % These specify the total degrees of freedom in each coordinate. 
            % It might be useful to seperate x/y expansions; given our current
            % assumptions, this is fine. 
            this.d_number_space   = number_cells_x(mesh_array{1});
            this.d_number_azimuth = number_azimuth(this.d_quadrature); 
            this.d_number_polar   = number_polar(this.d_quadrature);
                          
            % For now, just build the maximum basis.  Only for very fine meshes
            % would this get ridiculous.  We could also look at alternative
            % basis sets, like DCT's.
            this.d_basis_space    = DiscreteLP(this.d_number_space-1);
            this.d_basis_polar    = DiscreteLP(this.d_number_polar-1);
            this.d_basis_azimuth  = DiscreteLP(this.d_number_azimuth*2-1);         
            
            % Initialize solvers and related things for each mesh.  This way,
            % we don't need to reconstruct everything for each new keff.
            for i = 1:length(mesh_array)
                
                % Empty external source
                q_e = Source(this.d_mesh_array{i}, this.d_number_groups);
                
                % Create/initialize fission source.
                this.d_fission_source{i} = ...
                    FissionSource(this.d_state, ...
                                  this.d_mesh_array{i}, ...
                                  this.d_mat);
                initialize(this.d_fission_source{i});
                
                this.d_solver{i} = KrylovMG(this.d_input, ...
                                            this.d_state, ...
                                            this.d_boundary, ...
                                            this.d_mesh_array{i}, ...
                                            this.d_mat, ...
                                            this.d_quadrature, ...
                                            q_e, ...
                                            this.d_fission_source{i});
                              
                % Build fission and absorption vector                              
                build_vectors(this, this.d_mesh_array{i});
            end

        end
        

        % ======================================================================
        %> @brief Compure and return responses for new keff.
        %> @param   keff            Value of keff to interpolate responses.
        %> @return                  Responses    
        % ======================================================================    
        function [R, F, A, L] = get_responses(this, keff)   
            
            % Log initial time
            t = toc;
            
            % Return stored responses if keff was used one or two times ago.
%             if keff ~= this.d_keff_last
                 update(this, keff);
%             end
            
            R = this.d_R;
            F = this.d_F;
            A = this.d_A;
            L = this.d_L;
            
            % Keep this keff.
            this.d_keff_last = keff;
            
            % Add to timer
            this.d_time = this.d_time + (toc - t);
            
        end
        
        
        function this = update(this, keff)
            
            % Initial time
            t = toc;
            
            % Number of (non energy) expansion orders
            number_orders_space_angle = (this.d_order_space + 1) * ...
                                        (this.d_order_azimuth + 1) * ...
                                        (this.d_order_polar + 1);
                        
            % Degrees of freedom per surface                      
            number_orders_total = this.d_number_groups * ...
                                  (this.d_order_space + 1) * ...
                                  (this.d_order_azimuth + 1) * ...
                                  (this.d_order_polar + 1);
                              
            % Degrees of freedom per node, i.e. the size of its R-block.                  
            block_size = 4 * number_orders_total;
                              
                              
            % assumes 2-d, i.e. 4 faces
            this.d_R = zeros(block_size, block_size, this.d_number_nodes);
            this.d_F = zeros(block_size, this.d_number_nodes);
            this.d_A = zeros(block_size, this.d_number_nodes);
            this.d_L = zeros(block_size, 4, this.d_number_nodes);      
            
            number_responses = number_orders_total * this.d_number_nodes;
            count = 0;
            
            for i = 1:this.d_number_nodes 
                
                set_keff(this.d_solver{i},  keff);
                
                % Temporary storage of response
                this.d_coef    = cell(4, 1);
                this.d_coef{1} = zeros(number_orders_space_angle);
                this.d_coef{2} = zeros(number_orders_space_angle);
                this.d_coef{3} = zeros(number_orders_space_angle);
                this.d_coef{4} = zeros(number_orders_space_angle);
                
                this.d_fiss = zeros(number_orders_total, 1);
                this.d_abso = zeros(number_orders_total, 1);
                this.d_leak = zeros(number_orders_total, 4);
                    
                rf_index = 0;
                for g_o = 1:this.d_number_groups
                    for s_o = 0:this.d_order_space
                        for a_o = 0:this.d_order_azimuth
                            for p_o = 0:this.d_order_polar
                                
                                rf_index = rf_index + 1;
                                
                                % Setup and solve.  Reset fission after.
                                set_orders(this.d_bc, g_o, s_o, a_o, p_o);
                                reset(this.d_fission_source{i});
                                out = solve(this.d_solver{i});
                                
                                % Expand the coefficients
                                expand(this, rf_index, i);
                                count = count + 1;
                                time_remaining = (toc-t)/count * ...
                                    (number_responses-count);
                                if (get(this.d_input, 'rf_print_out'))
                                    fprintf(' ResponseDriver: i=%2i g_o=%2i s_o=%2i a_o=%2i p_o=%2i trem=%12.8f\n', ...
                                        i, g_o, s_o, a_o, p_o, time_remaining );
                                end
                                
                            end % azimuth loop
                        end % polar loop
                    end % space loop
                end % group loop

                % Build the full response matrix.  In the future, this can be
                % generalized to have incident conditions on all surfaces.
                mo = number_orders_total;
                R( (0*mo)+1: 1*mo, (0*mo)+1: 1*mo) = this.d_coef{1}(:, :); % left -> left
                R( (1*mo)+1: 2*mo, (0*mo)+1: 1*mo) = this.d_coef{2}(:, :); % left -> right
                R( (2*mo)+1: 3*mo, (0*mo)+1: 1*mo) = this.d_coef{3}(:, :); % left -> top
                R( (3*mo)+1: 4*mo, (0*mo)+1: 1*mo) = this.d_coef{4}(:, :); % left -> bottom
                
                R( (0*mo)+1: 1*mo, (1*mo)+1: 2*mo) = this.d_coef{2}(:, :); % right -> lef
                R( (1*mo)+1: 2*mo, (1*mo)+1: 2*mo) = this.d_coef{1}(:, :); % right -> right
                R( (2*mo)+1: 3*mo, (1*mo)+1: 2*mo) = this.d_coef{4}(:, :); % right -> top
                R( (3*mo)+1: 4*mo, (1*mo)+1: 2*mo) = this.d_coef{3}(:, :); % right -> bottom
                
                R( (0*mo)+1: 1*mo, (2*mo)+1: 3*mo) = this.d_coef{4}(:, :); % bottom -> left
                R( (1*mo)+1: 2*mo, (2*mo)+1: 3*mo) = this.d_coef{3}(:, :); % bottom -> right
                R( (2*mo)+1: 3*mo, (2*mo)+1: 3*mo) = this.d_coef{1}(:, :); % bottom -> bottom
                R( (3*mo)+1: 4*mo, (2*mo)+1: 3*mo) = this.d_coef{2}(:, :); % bottom -> top
                
                R( (0*mo)+1: 1*mo, (3*mo)+1: 4*mo) = this.d_coef{3}(:, :); % top -> left
                R( (1*mo)+1: 2*mo, (3*mo)+1: 4*mo) = this.d_coef{4}(:, :); % top -> right
                R( (2*mo)+1: 3*mo, (3*mo)+1: 4*mo) = this.d_coef{2}(:, :); % top -> bottom
                R( (3*mo)+1: 4*mo, (3*mo)+1: 4*mo) = this.d_coef{1}(:, :); % top -> top
                
                L1 = this.d_leak(:, 1)'; % into 1 ->  leak from 1
                L2 = this.d_leak(:, 2)'; % into 1 ->  leak from 2
                L3 = this.d_leak(:, 3)'; % etc. Repeat for assumed sym.
                L4 = this.d_leak(:, 4)';
                
                this.d_L(:, 1, i) = [L1 L2 L4 L3]'; %   |J1|   leak from side 1
                this.d_L(:, 2, i) = [L2 L1 L3 L4]'; % * |J2| = leak from side 2
                this.d_L(:, 3, i) = [L3 L4 L1 L2]'; %   |J3|   etc.
                this.d_L(:, 4, i) = [L4 L3 L2 L1]'; %   |J4|
                
                this.d_R(:, :, i) = R;
                this.d_F(:,    i) = [this.d_fiss; this.d_fiss; this.d_fiss; this.d_fiss];
                this.d_A(:,    i) = [this.d_abso; this.d_abso; this.d_abso; this.d_abso;];
                
            end % node loop
            t_final = toc-t;
            fprintf(' ResponseDriver total time (seconds): %12.8f\n', t_final);
        end
        
        function this = expand(this, index_in, node_idx)
            
            q_f = this.d_fission_source{node_idx};
            num_ang  = this.d_number_polar*this.d_number_azimuth;
            
            % Octants for incident left conditions
            octants = [ 3 2   % ref
                        1 4   % far
                        4 3   % lef
                        2 1]; % rig
            
            % Expand the coefficients
            for side = 1:4
                
                f = zeros(this.d_number_space,num_ang,this.d_number_groups);
              
                for g = 1:this.d_number_groups
                    
                    for o = 1:2
                        o_in  = octants(side, o); % incident octant
                        % always left to right in angle w/r to outgoing flux
                        if o == 1
                            a1 = 1;
                            a2 = this.d_number_azimuth;
                            a3 = 1;
                        else
                            a1 = this.d_number_azimuth;
                            a2 = 1;
                            a3 = -1;
                        end
                        % Get psi(x, angles)
                        set_group(this.d_boundary, g);
                        if (side == 1 || side == 2)
                            ft = get_psi_v_octant(this.d_boundary, o_in, Boundary.OUT);
                            if (side == 1)
                                ft(1:end,:)=ft(end:-1:1,:); % reverse space
                            end
                        else
                            ft = get_psi_h_octant(this.d_boundary, o_in, Boundary.OUT);
                            if (side == 4)
                                ft(1:end,:)=ft(end:-1:1,:);  % reverse space
                            end
                        end
                        % populate the vectors we expand
                        ang = 1;
                        for a = a1:a3:a2
                            for p = 1:this.d_number_polar
                                f(:, (o-1)*num_ang+ang,g) =  ...
                                    ft(:, (a-1)*this.d_number_polar+p);
                                ang = ang + 1;
                            end
                        end
                    end
                    
                end
                
                % Group->Space->Azimuth->Polar.  f(space, angle, group)
                index_out = 1;
                for g = 1:this.d_number_groups
                    for ord_s = 1:this.d_order_space+1
                        for ord_a = 1:this.d_order_azimuth+1
                            
                            
                            b = this.d_basis_azimuth(:, ord_a);
                            % Reorder if it's a horizontal surface, since
                            % [mu,eta] are always in that order
                            if side > 2
                                lb = length(b);
                                b(1:lb/2) = b(lb/2:-1:1);
                                b(lb/2+1:end) = b(end:-1:lb/2+1);
                            end
                            
                            for ord_p = 1:this.d_order_polar+1
                                tmp = f(:, :, g)'*this.d_basis_space(:, ord_s);                        
                                psi_ap = reshape(tmp,this.d_number_polar,2*this.d_number_azimuth);
                                psi_p = psi_ap*b;
                                this.d_coef{side}(index_out, index_in) = ...
                                    psi_p'*this.d_basis_polar(:, ord_p); % i <- k
                                index_out = index_out + 1;
                            end % out p
                            
                        end % out a
                    end % out s
                end % out g
                
            end % side loop
            
            % Fission and Absorption rates
            update(q_f);
           
            for g = 1:number_groups(this.d_mat)
                phi = flux(this.d_state, g);
                this.d_fiss(index_in) = ...
                    this.d_fiss(index_in) + ...
                    sum(this.d_volume{node_idx}.*phi.*this.d_fission{node_idx}(:, g));
                this.d_abso(index_in) = ...
                    this.d_abso(index_in) + ...
                    sum(this.d_volume{node_idx}.*phi.*this.d_absorption{node_idx}(:, g));    
                set_group(this.d_boundary, g);
                this.d_leak(index_in,:) = ... 
                    this.d_leak(index_in,:) + ...
                    get_leakage(this.d_boundary);
            end

            % Pin powers
            % to be added
            
        end
        
    end
    
    methods (Access = private)
        
        function this = build_vectors(this, mesh) 
            
            idx = length(this.d_fission)+1;
            
            fission    = zeros(number_cells(mesh), number_groups(this.d_mat));
            absorption = zeros(number_cells(mesh), number_groups(this.d_mat));
            volume     = zeros(number_cells(mesh), 1);
            
            % Get the fine mesh material map or the region material map.
            n  = number_cells(mesh);
            ng = number_groups(this.d_mat);
            if meshed(mesh)
                mat = reshape(mesh_map(mesh, 'MATERIAL'), n, 1);
            else                
                mat = region_mat_map(mesh);
            end            
            for i = 1:n
                for g = 1:ng
                    fission(i, g)    = nu_sigma_f(this.d_mat, mat(i), g);
                    absorption(i, g) = sigma_a(this.d_mat, mat(i), g);
                end
            end
            k = 1;
            for j = 1:number_cells_y(mesh)
                for i = 1:number_cells_x(mesh)
                    volume(k) = dx(mesh, i)*dy(mesh, j);
                    k = k + 1;
                end
            end
            
            this.d_fission{idx} = fission;
            this.d_absorption{idx} = absorption;
            this.d_volume{idx} = volume;

        end
        
        
    end
    
end