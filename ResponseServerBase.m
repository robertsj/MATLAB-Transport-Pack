%> @file  ResponseServerBase.m
%> @brief ResponseServerBase class definition.
% ==============================================================================
%> @brief Base class for serving response functions.
%
%> Relevant input database entries include
%>   - rf_order_space
%>   - rf_order_azimuth
%>   - rf_order_polar
%>   - dimension
%> which specify the orders used in the computation.
%>
%> @todo I've attempted no handling of order reduction, e.g. 2nd order means
%>       being limited to total orders <= 2, e.g. 
%>       order(space)*order(angle) <= 2.
% ==============================================================================
classdef ResponseServerBase < handle

    properties (Access = public)
        %> User input database.
        d_input
        %> Time spent in the server.
        d_time
        %> Problem dimension
        d_dimension
        %> Number of nodes
        d_number_nodes
        %> Number of energy groups
        d_number_groups           
        %> Working spatial order
        d_order_space
        %> Working azimuthal order
        d_order_azimuth
        %> Working polar order
        d_order_polar
        %> Maximum spatial order
        d_max_order_space
        %> Maximum azimuthal order
        d_max_order_azimuth
        %> Maximum polar order
        d_max_order_polar        
        %> Response matrix from most recent update
        d_R
        %> Fission operator from most recent update
        d_F
        %> Absorption operator from most recent update
        d_A
        %> Leakage operator from most recent update
        d_L        
    end
    
    methods
        
        % ======================================================================
        %> @brief  Class constructor.
        %> @param  input    Input database.
        %> @return          Instance of the ResponseServerBase class.
        % ======================================================================
        function this = ResponseServerBase(input) 
            this.d_input = input;
            
            % Read in the basic information
            this.d_dimension        = get(this.d_input, 'dimension');
            this.d_number_groups    = get(this.d_input, 'number_groups');
            this.d_order_space      = get(this.d_input, 'rf_order_space');
            this.d_order_azimuth    = get(this.d_input, 'rf_order_azimuth');
            this.d_order_polar      = get(this.d_input, 'rf_order_polar');    
            this.d_number_nodes     = get(this.d_input, 'rf_number_nodes');
            
            % By default, the maximum order *is* the working order.  In cases
            % where operators are already computed to high order, the client
            % can set a low order and get the corresponding operators.
            this.d_max_order_space   = this.d_order_space;
            this.d_max_order_azimuth = this.d_order_azimuth;
            this.d_max_order_polar   = this.d_order_polar;
        end
        
        % ======================================================================
        %> @brief Return responses for the given eigenvalue.
        %
        %> The arrays R, F, A, and L have dimensions @f$(M,M,N)@f$,  
        %> @f$(M,N)@f$,  @f$(M,N)@f$,  and @f$(M,S,N)@f$,  respectively, where
        %> @f$ M @f$, is the number of degrees of freedom per node,
        %> @f$ N @f$ is the number of unique nodes, and @f$ S @f$ is the 
        %> number of surfaces per node.
        %>
        %> @param   keff 	Global eigenvalue defining local problems.
        %> @return          Responses    
        % ======================================================================    
        [R, F, A, L] = get_responses(this, keff)   

        
        %> @name Getters
        %> @{

        function n = number_nodes(this)
            n = this.d_number_nodes;
        end    

        function t = time(this)
            t = this.d_time;
        end
        
        %> @}
        
        %> @name Setters
        %> @{
        
        function set_orders(this, orders)
           
            if isfield(orders, 'mso')
                this.d_max_order_space = orders.mso;
            end
            if isfield(orders, 'so')
                this.d_order_space = orders.so;
            end
            if isfield(orders, 'mao')
                this.d_max_order_azimuth = orders.mao;
            end
            if isfield(orders, 'ao')
                this.d_order_azimuth = orders.ao;
            end
            if isfield(orders, 'mpo')
                this.d_max_order_polar = orders.mpo;
            end
            if isfield(orders, 'po')
                this.d_order_polar = orders.po;
            end            
            
        end
        
        %> @}        
        
    end
    
    methods (Access = protected)
                
        % ======================================================================
        %> @brief  Extract low order operators from higher order operators.
        %
        %> The low orders are based on the "working" orders whereas the actual
        %> (high) orders are based on previous computation or database.
        %>
        %> @param  R        High order response matrix.
        %> @param  F        High order fission operator.
        %> @param  A        High order absorption operator.
        %> @param  L        High order leakage operator.
        %> @return          Low order operators.
        % ======================================================================
        function [R2, F2, A2, L2] = reduce(this, R, F, A, L)
            
            % Go through and build an index of all order permutations that meet
            % the requested order combination.
            n = 0;
            row   = 1;
            for side = 1:2*this.d_dimension
                for g = 1:this.d_number_groups
                    for s = 0:this.d_max_order_space
                        for a = 0:this.d_max_order_azimuth
                            for p = 0:this.d_max_order_polar
                                if (s <= this.d_order_space   && ...
                                    a <= this.d_order_azimuth && ...
                                    p <= this.d_order_polar)
                                    n = n + 1;
                                    idx(n) = row;
                                end
                                row = row + 1;
                            end
                        end
                    end
                end
            end
            R2(1:n,1:n,:,:) = R(idx, idx, :, :);
            L2(1:n,:,:,:)   = L(idx, :, :, :);
            F2(1:n,:,:)     = F(idx, :, :);
            A2(1:n,:,:)     = A(idx, :, :);
        end
        
    end
        
    
end