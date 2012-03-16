%> @file  ResponseDriver.m
%> @brief Drives generation of response functions.
% ==============================================================================
%> @brief ResponseDriver class definition.
%
%> This drives generation of response functions for a given set of local
%> problems (aka "elements" or "nodes").  The user specifies a set of values
%> for keff and the maximum orders in space and angle.  The output is placed
%> into an HDF5 file.
%> 
%> In this initial implementation, we assume that all nodes have the same
%> meshing, and that all nodes are symmetric.  Hence, only the LEFT faces is
%> used for incident conditions.
%> 
% ==============================================================================
classdef ResponseDriver < handle

    properties (Access = private)
        d_input
        d_mat
        d_mesh_array
        d_number_nodes
        d_solver_array
        d_fission_source
        d_quadrature
        d_boundary
        d_bc
        d_state
        d_R
        d_F
        d_A
        d_L
        d_fiss
        d_abso
        d_leak
        d_fission
        d_absorption
        d_volume
        d_coef
        d_max_o
        d_max_g_o
        d_max_s_o
        d_max_a_o
        d_max_p_o
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
        d_k_vector
    end
    
    properties (Access = public)
       rates 
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function this = ResponseDriver(input, mat, mesh_array)
            
            disp('*** Building ResponseDriver')
            
            % Verify input and materials have same group count.
            assert(get(input, 'number_groups') == number_groups(mat));
            % Assert we got at least one mesh
            assert(~isempty(mesh_array));
            
            % Force the boundaries to be be what we can handle.
            put(input, 'bc_left',               'response');
            put(input, 'bc_right',              'vacuum');
            put(input, 'bc_top',                'vacuum');
            put(input, 'bc_bottom',             'vacuum');
            
            this.d_input = input;
            this.d_mat   = mat;
            this.d_mesh_array = mesh_array;
            
            % Setup the quadrature.  Could default to 2, but be explicit!
            qo = get(input, 'quad_order');
            if (~qo)
                error('user must specify a QuadrupleRange order')
            end
            this.d_quadrature  = QuadrupleRange(qo);
            
            % I'm not sure how to generalize for MOC at the moment.
            this.d_boundary = ...
                BoundaryMesh(input, mesh_array{1}, this.d_quadrature);
            % Hackish assumption: only incident on left.
            this.d_bc       = get_bc(this.d_boundary, Mesh.LEFT);
            

            this.d_number_nodes = length(mesh_array);
            this.d_state  = State(input, mesh_array{1});
            
            this.d_number_space   = number_cells_x(mesh_array{1});
            this.d_number_azimuth = number_azimuth(this.d_quadrature); 
            this.d_number_polar   = number_polar(this.d_quadrature);
                          
            this.d_basis_space    = DiscreteLP(this.d_number_space-1);
            this.d_basis_polar    = DiscreteLP(this.d_number_polar-1);
            this.d_basis_azimuth  = DiscreteLP(this.d_number_azimuth*2-1);            

        end
        
        % ======================================================================
        %> @brief Run
        % ======================================================================
        function this = run(this)
            tic
            disp('running...')
            this.d_max_g_o = get(this.d_input, 'number_groups');
            this.d_max_s_o = get(this.d_input, 'rf_max_order_space');
            this.d_max_a_o = get(this.d_input, 'rf_max_order_azimuth');
            this.d_max_p_o = get(this.d_input, 'rf_max_order_polar');
         
            this.d_k_vector = get(this.d_input, 'rf_k_vector');
            if (~this.d_k_vector)
                this.d_k_vector = 1.0;
            end
            
            % Empty external source
            q_e = Source(this.d_mesh_array{1}, this.d_max_g_o);
            
            this.d_max_o = (this.d_max_s_o + 1) * ...
                 (this.d_max_a_o + 1)*(this.d_max_p_o + 1);
            mo = this.d_max_g_o * (this.d_max_s_o + 1) * ...
                 (this.d_max_a_o + 1)*(this.d_max_p_o + 1);     
             
            % assumes 2-d, i.e. 4 faces
            this.d_R = zeros(4*mo, 4*mo, length(this.d_k_vector), this.d_number_nodes);
            this.d_F = zeros(4*mo,       length(this.d_k_vector), this.d_number_nodes);
            this.d_A = zeros(4*mo,       length(this.d_k_vector), this.d_number_nodes);
            this.d_L = zeros(4*mo, 4,    length(this.d_k_vector), this.d_number_nodes);      
            
            number_responses = mo * this.d_number_nodes * ...
                               length(this.d_k_vector);
            count = 0;
            for i = 1:this.d_number_nodes
                
                q_f = FissionSource(this.d_state, ...
                                    this.d_mesh_array{i}, ...
                                    this.d_mat);
                initialize(q_f);
                
                solver = KrylovMG(this.d_input, ...
                                  this.d_state, ...
                                  this.d_boundary, ...
                                  this.d_mesh_array{i}, ...
                                  this.d_mat, ...
                                  this.d_quadrature, ...
                                  q_e, ...
                                  q_f);
                              
                % Build fission and absorption vector                              
                build_vectors(this, this.d_mesh_array{i});
                             
                for ki = 1:length(this.d_k_vector)
                
                    set_keff(solver,  this.d_k_vector(ki)); 

                    % Temporary storage of response
                    this.d_coef    = cell(4, 1);
                    this.d_coef{1} = zeros(this.d_max_o);
                    this.d_coef{2} = zeros(this.d_max_o);
                    this.d_coef{3} = zeros(this.d_max_o);
                    this.d_coef{4} = zeros(this.d_max_o);
                    
                    this.d_fiss = zeros(mo, 1);
                    this.d_abso = zeros(mo, 1);
                    this.d_leak = zeros(mo, 4);
                    
                    rf_index = 0;
                    for g_o = 1:this.d_max_g_o
                        for s_o = 0:this.d_max_s_o
                            for a_o = 0:this.d_max_a_o
                                for p_o = 0:this.d_max_p_o

                                    rf_index = rf_index + 1;

                                    % Setup and solve.  Reset fission after.
                                    set_orders(this.d_bc, g_o, s_o, a_o, p_o);
                                    reset(q_f);
                                    out = solve(solver);
                                    
                                    % Expand the coefficients
                                    expand(this, rf_index, q_f);
%                                     rates = reaction_rates(this.d_mat, ... 
%                                         this.d_mesh_array{i}, this.d_state, ...
%                                         this.d_boundary)
%                                     this.rates = rates;
                                    count = count + 1;
                                    time_remaining = toc/count * ...
                                        (number_responses-count);
                                    fprintf(' ResponseDriver: i=%2i ki=%2i g_o=%2i s_o=%2i a_o=%2i p_o=%2i trem=%12.8f\n', ...
                                        i, ki, g_o, s_o, a_o, p_o, time_remaining );

                                end % azimuth loop
                            end % polar loop
                        end % space loop
                    end % group loop

                    % Build the full response matrix.  In the future, this can be
                    % generalized to have incident conditions on all surfaces.

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
         
                    this.d_L(:, 1, ki, i) = [L1 L2 L4 L3]'; %   |J1|   leak from side 1
                    this.d_L(:, 2, ki, i) = [L2 L1 L3 L4]'; % * |J2| = leak from side 2
                    this.d_L(:, 3, ki, i) = [L3 L4 L1 L2]'; %   |J3|   etc.
                    this.d_L(:, 4, ki, i) = [L4 L3 L2 L1]'; %   |J4|

                    this.d_R(:, :, ki, i) = R;
                    this.d_F(:,    ki, i) = [this.d_fiss; this.d_fiss; this.d_fiss; this.d_fiss];
                    this.d_A(:,    ki, i) = [this.d_abso; this.d_abso; this.d_abso; this.d_abso;];
                end % kloop
                
            end
            t_final = toc;
            fprintf(' ResponseDriver total time (seconds): %12.8f\n', t_final);
        end

        % ======================================================================
        %> @brief Class constructor
        % ======================================================================        
        function [R, F, A, L] = update(this, keff)
            
        end
        
        function this = expand(this, index_in, q_f)
            
            
            num_ang        = this.d_number_polar*this.d_number_azimuth;
            
            octants = [ 3 2   % ref
                        1 4   % far
                        4 3   % lef
                        2 1]; % rig
            
            % Expand the coefficients
            for side = 1:4
                
                f = zeros(this.d_number_space,num_ang,this.d_max_g_o);
                
                for g = 1:this.d_max_g_o
                    
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
                        %for s = 1:this.d_number_space
                            ang = 1;
                            for a = a1:a3:a2
                                for p = 1:this.d_number_polar
                                    f(:, (o-1)*num_ang+ang,g) =  ...
                                        ft(:, (a-1)*this.d_number_polar+p);
                                    ang = ang + 1;
                                end
                            end
                        %end
                    end
                    
                end
                
                % Group->Space->Azimuth->Polar.  f(space, angle, group)
                
                
                index_out = 1;
                for g = 1:this.d_max_g_o
                    for ord_s = 1:this.d_max_s_o+1
                        for ord_a = 1:this.d_max_a_o+1
                            
                            b = this.d_basis_azimuth(:, ord_a);
                            if side > 2
                                lb = length(b);
                                b(1:lb/2) = b(lb/2:-1:1);
                                b(lb/2+1:end) = b(end:-1:lb/2+1);
                            end
                            
                            for ord_p = 1:this.d_max_p_o+1
                               % psi_ap = zeros(2*this.d_number_azimuth, this.d_number_polar);
                               % angle = 0;
                                tmp = f(:, :, g)'*this.d_basis_space(:, ord_s);
                                %az=0;
%                                 for o = 1:2
%                                     a1 = 1;
%                                     a2 = this.d_number_azimuth;
%                                     a3 = 1;
%                                     for a = a1:a3:a2
%                                         az = az+1;
%                                         for p = 1:this.d_number_polar
%                                             angle = angle + 1;
%                                             psi_ap(az, p) = tmp(angle);%f(:, angle, g)'*this.d_basis_space(:, ord_s);
%                                         end
%                                     end
%                                 end
                                
%                                 tmp2 = reshape(tmp(1:end/2),this.d_number_polar,[]);
%                                 psi_ap(1:this.d_number_azimuth,1:this.d_number_polar) = tmp2';
%                                 tmp2 = reshape(tmp(end/2+1:end),this.d_number_polar,[]);
%                                 psi_ap(this.d_number_azimuth+1:end,1:this.d_number_polar) = tmp2';
%                                 
                               % psi_ap2 = psi_ap;
                                psi_ap = reshape(tmp,this.d_number_polar,2*this.d_number_azimuth);
                                %psi_ap = tmp2';
                                
                                
                                
                                
%                                 psi_p = zeros(this.d_number_polar, 1);
%                                 for p = 1:this.d_number_polar
%                                     psi_p(p) = psi_ap(:, p)*b; % [np,1]=[np,na][na,1]
%                                 end
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
                this.d_fiss(index_in) = this.d_fiss(index_in) + ...
                    sum(this.d_volume.*phi.*this.d_fission(:, g));
                this.d_abso(index_in) = this.d_abso(index_in) + ...
                    sum(this.d_volume.*phi.*this.d_absorption(:, g));    
                
                set_group(this.d_boundary, g);
                this.d_leak(index_in,:) = this.d_leak(index_in,:) + ...
                    get_leakage(this.d_boundary);
            end

            % Pin powers
            % to be added
            
        end
        
        function [R, F, A, L] = get_responses(this)
           R = this.d_R;
           F = this.d_F;
           A = this.d_A;
           L = this.d_L;
        end
        
    end
    
    methods (Access = private)
        
        function this = build_vectors(this, mesh) 
            this.d_fission    = zeros(number_cells(mesh), number_groups(this.d_mat));
            this.d_absorption = zeros(number_cells(mesh), number_groups(this.d_mat));
            this.d_volume     = zeros(number_cells(mesh), 1);
            % Get the fine mesh material map or the region material map.
            n = number_cells(mesh);
            ng = number_groups(this.d_mat);
            if meshed(mesh)
                mat = reshape(mesh_map(mesh, 'MATERIAL'), n, 1);
            else                
                mat = region_mat_map(mesh);
            end            
            for i = 1:n
                for g = 1:ng
                    this.d_fission(i, g)    = nu_sigma_f(this.d_mat, mat(i), g);
                    this.d_absorption(i, g) = sigma_a(this.d_mat, mat(i), g);
                end
            end
            k = 1;
            for i = 1:number_cells_x(mesh)
                for j = 1:number_cells_y(mesh)
                    this.d_volume(k) = dx(mesh, i)*dy(mesh, j);
                    k = k + 1;
                end
            end
        end
        
    end
    
end