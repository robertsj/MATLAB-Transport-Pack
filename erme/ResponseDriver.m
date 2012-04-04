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
            this.d_max_order_space   = input.get('rf_max_order_space');
            this.d_max_order_azimuth = input.get('rf_max_order_azimuth');
            this.d_max_order_polar   = input.get('rf_max_order_polar'); 
            this.d_keff_vector       = input.get('rf_keff_vector');
            this.d_node_descriptions = input.get('rf_node_descriptions');
            if length(this.d_node_descriptions) ~= this.d_number_nodes
                warning('user:input', ...
                    'Wrong number or missing node descroptions; using default for all.')
                this.d_node_descriptions = cell(this.d_number_nodes, 1);
                for i = 1:this.d_number_nodes
                    this.d_node_descriptions{:} = ['assembly', num2str(i)];
                end
            end
            
            % Create the server
            this.d_server = ResponseServer(input, mat, mesh_array);
            % Set the requested orders from the server to be the maximum the
            % database shall hold.
            orders.so = this.d_max_order_space;
            orders.ao = this.d_max_order_azimuth;
            orders.po = this.d_max_order_polar;
            set_orders(this.d_server, orders);
            
            % Display
            tot = (this.d_max_order_space+1)*(this.d_max_order_azimuth+1)*...
                   (this.d_max_order_polar+1)*(this.d_number_groups);
            fprintf(' RESPONSE DRIVER \n'); %, Error: %12.8f
            fprintf('     number of nodes: %5i \n', length(mesh_array));
            fprintf('     number of keffs: %5i \n', length(this.d_keff_vector));
            fprintf('    number of groups: %5i \n', this.d_number_groups);
            fprintf('       spatial order: %5i \n', this.d_max_order_space);
            fprintf('     azimuthal order: %5i \n', this.d_max_order_azimuth );
            fprintf('         polar order: %5i \n', this.d_max_order_polar);
            fprintf('*** number responses: %5i \n', tot);
        end
        
        % ======================================================================
        %> @brief Run
        % ======================================================================
        function this = run(this)
            
            t = toc;
            
            % Total number of nodes * keffs
            total = 0;
            for n = 1:this.d_number_nodes
            	total = total + length(this.d_keff_vector{n});
            end
            
            count = 0;
            
            for node_index = 1:this.d_number_nodes
            
                for k_index = 1:length(this.d_keff_vector{node_index})
                 
                    keff = this.d_keff_vector{node_index}(k_index);
                    
                    
                    [R, F, A, L] = ...
                        get_node_responses(this.d_server, node_index, keff); 
                    
                    this.d_R_all{node_index}(:, :, k_index) = R(:, :);
                    this.d_L_all{node_index}(:, :, k_index) = L(:, :);
                    this.d_F_all{node_index}(:,    k_index) = F(:);
                    this.d_A_all{node_index}(:,    k_index) = A(:);
                    
                    count = count + 1;
                    time_remaining = (toc-t)/count*(total-count);
                    
                    if (get(this.d_input, 'rf_print_out'))
                        fprintf(' ResponseDriver.run: n=%2i k=%2i trem=%12.8f\n', ...
                            node_index, k_index, time_remaining );
                    end
                    
                end
                
            end
            
            this.initialize_write();
            
            for node_index = 1:this.d_number_nodes
                write_response(this, ...
                               node_index, ...
                               this.d_node_descriptions{node_index}, ...
                               this.d_R_all{node_index}(:,:,:), ...
                               this.d_F_all{node_index}(:,:), ...
                               this.d_A_all{node_index}(:,:), ...
                               this.d_L_all{node_index}(:,:,:));
            end
            
            this.d_time = this.d_time + (toc - t);
            if (get(this.d_input, 'rf_print_out'))
                fprintf(' ResponseDriver.run total time (seconds): %12.8f\n', this.d_time);
            end
        end

        
    end
    
   
end