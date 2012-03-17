%> @file  ResponseDriver.m
%> @brief ResponseDriver class driver class definition.
% ==============================================================================
%> @brief Drives generation of a response function database.
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
classdef ResponseDriver < ResponseDB

    properties (Access = private)      
        %> Response function server (i.e. the generator)
        d_server
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function this = ResponseDriver(input, mat, mesh_array)
            
            % Call base class
            this = this@ResponseDB(input);
            DBC.Require('this.number_nodes()==length(mesh_array)');
            this.d_max_order_space   = get(this.d_input, 'rf_max_order_space');
            this.d_max_order_azimuth = get(this.d_input, 'rf_max_order_azimuth');
            this.d_max_order_polar   = get(this.d_input, 'rf_max_order_polar'); 
            this.d_keff_vector = input.get('rf_keff_vector');
            
            % Create the server
            this.d_server = ResponseServer(input, mat, mesh_array);
            % Set the requested orders from the server to be the maximum the
            % database shall hold.
            orders.so = this.d_max_order_space;
            orders.ao = this.d_max_order_azimuth;
            orders.po = this.d_max_order_polar;
            set_orders(this.d_server, orders);
            
            
        end
        
        % ======================================================================
        %> @brief Run
        % ======================================================================
        function this = run(this)
            
            t = toc;
 
            for ki = 1:length(this.d_keff_vector)
                
                [R, F, A, L] = get_responses(this.d_server, ...
                                             this.d_keff_vector(ki)); 
                
                for i = 1:this.d_number_nodes
                    this.d_R_all(:, :, ki, i) = R(:, :, i);
                    this.d_L_all(:, :, ki, i) = L(:, :, i);
                    this.d_F_all(:,    ki, i) = F(:, i);
                    this.d_A_all(:,    ki, i) = A(:, i);
                end
                
            end
            
            this.initialize_write();
            
            for i = 1:this.d_number_nodes
                write_response(this, i, ['assembly',num2str(i)], ...
                               this.d_R_all(:,:,:,i), ...
                               this.d_F_all(:,:,i), ...
                               this.d_A_all(:,:,i), ...
                               this.d_L_all(:,:,:,i));
            end
            
            this.d_time = this.d_time + (toc - t);
            if (get(this.d_input, 'rf_print_out'))
                fprintf(' ResponseDriver total time (seconds): %12.8f\n', this.d_time);
            end
        end

        
    end
    
   
end