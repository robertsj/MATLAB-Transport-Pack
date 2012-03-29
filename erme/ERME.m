%> @file  ERME.m
%> @brief ERME class definition.
% ==============================================================================
%> @brief Solve a problem via the eigenvalue response matrix method.
% ==============================================================================
classdef ERME < handle
    
    properties (Access = private)
        
        d_input
        
        %> @name Geometry
        %> @{        
        
        %> Connectivity matrix.
        d_connect
        d_M
        d_leak
        
        %> Stencil dimensions
        d_dimension = 0;
        d_number_x = 0;
        d_number_y = 0;
        d_number_z = 0;
        
        %> List of element types
        d_elements;
        
        %> @}
        
        %> @name Responses
        %> @{
        
        %> Response matrix.
        d_R
        %> Fission operator.
        d_F
        %> Absorption operator.
        d_A
        %> Leakage operator.
        d_L
        
        %> Size of unknown
        d_size
        
        %> Number of elements
        d_number_elements
        
        %> Response function server
        d_server
        
        %> @}
        
        %> @name State
        %> @{
        
        %> Incident current coefficients.
        d_J
        %> k-Eigenvalue.
        d_keff
        %> lambda-Eigenvalue.
        d_lambda
        
        %> @}
        
    end
    
    methods
        
        % ======================================================================
        %> @brief  Class constructor
        %> @param  input    Input database.
        %> @param  server   Response function server.
        %> @param  geometry Problem geo
        %> @return Instance of the ERME class.
        % ======================================================================
        function this = ERME(input, server, elements)
            this.d_input = input;
            % The server is either
            this.d_server = server;
            % Check requested orders against the server's max
%             if (get(input, 'rf_space_order') > max_space_order(server));
%                 put(input, 'rf_space_order',  max_space_order(server));
%             end
%             if (get(input, 'rf_azimuth_order') > max_azimuth_order(server));
%                 put(input, 'rf_azimuth_order',  max_azimuth_order(server));
%             end
%             if (get(input, 'rf_polar_order') > max_polar_order(server));
%                 put(input, 'rf_polar_order',  max_polar_order(server));
%             end            
            % The elements are arranged "physically" as input. Transpose.
            for k = 1:length(elements(1, 1, :))
                tmp(:, :) = elements(:, :, k);
                elem(:, :, k) = flipud(tmp)';
            end
            elements = elem;
            %
            this.d_dimension = get(input, 'dimension');

            % Stencil size and element count.
            this.d_number_x = length(elements(:, 1, 1));
            this.d_number_y = length(elements(1, :, 1));
            this.d_number_z = length(elements(1, 1, :));
            this.d_number_elements = sum(sum(sum(elements>0)));
            
            % size of current = (number faces) * 
                 
                      
            % Generate the connectivity
            this.d_connect = Connect(input, elements);  
            [this.d_M this.d_leak] = build(this.d_connect);
            
            % Size of J
            this.d_size = length(this.d_M);
            
            
            % Define element type vector
            this.d_elements = zeros(this.d_number_elements, 1);
            e_idx = 0;
            for k = 1:this.d_number_z
                for j = 1:this.d_number_y
                    for i = 1:this.d_number_x                
                        if elements(i, j, k) > 0
                            e_idx = e_idx + 1;
                            this.d_elements(e_idx) = elements(i, j, k);
                        end
                        
                    end
                end
            end
            assert(e_idx == this.d_number_elements);
        end
        
        % ======================================================================
        %> @brief  Update the responses for a new eigenvalue.
        %> @param  keff     Updated eigenvalue.
        % ======================================================================        
        function update(this, keff)
            % Get new data.
            [R, F, A, L] = get_responses(this.d_server, keff);
           
            % Initialize.
            R_list = cell(this.d_number_elements, 1);
            L_list = cell(this.d_number_elements, 1);
            this.d_F = zeros(1, this.d_size);
            this.d_A = zeros(1, this.d_size);
            
            % Assemble.
            size_block = length(R(:, 1, 1, 1));
            for i = 1:this.d_number_elements
                type  = this.d_elements(i);
                range = ((i-1)*size_block+1):(i*size_block);
                tmp1(1:size_block, 1:size_block) = R(:, :, type);
                R_list{i} = sparse(tmp1);
                this.d_F(1, range) = F(:, type);
                this.d_A(1, range) = A(:, type);
                tmp2(1:size_block, 1:4) = L(:, 1:4, type);
                L_list{i} = sparse(tmp2');
            end
            this.d_R = blkdiag(R_list{:});
            this.d_L = blkdiag(L_list{:});
        end
        
        % ======================================================================
        %> @brief  Update the responses for a new eigenvalue.
        %> @param  keff     Updated eigenvalue.
        % ======================================================================        
        function update_state(this, J, keff, lambda)
            this.d_J      = J;
            this.d_keff   = keff;
            this.d_lambda = lambda;
        end
        
        % ======================================================================
        %> @brief  Get the operators.
        % ======================================================================        
        function [R, F, A, L, M, leak] = get_operators(this)
            R = this.d_R;
            F = this.d_F;
            A = this.d_A;
            L = this.d_L;
            M = this.d_M;
            leak = this.d_leak;
        end
        
        % ======================================================================
        %> @brief  Return connectivity.
        % ======================================================================    
        function M = get_connect(this)
            M = this.d_connect; 
        end
        
        % ======================================================================
        %> @brief  Return connectivity.
        % ======================================================================    
        function n = number_elements(this)
            n = this.d_number_elements; 
        end        
        
        % ======================================================================
        %> @brief  Return surfaces per elements.
        % ======================================================================    
        function n = number_faces(this)
            n = 2*this.d_dimension; 
        end                
        
        % ======================================================================
        %> @brief  Return size of unknown.
        % ======================================================================    
        function n = size_J(this)
            n = this.d_size; 
        end
        
        function k = get_keff(this)
            k = this.d_keff;
        end
        function l = get_lambda(this)
            l = this.d_lambda;
        end
        
    end
    
end