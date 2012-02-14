%> @file  FLARE.m
%> @brief FLARE class definition.
% ==============================================================================
%> @brief  Implements the classic FLARE method in 2-D.
%
%> More here.
% ==============================================================================
classdef FLARE < handle
    
    properties
        %> User input.
        d_input                 
        %> Material database.
        d_mat                   
        %> Core map (2-d array of bundle id's)
        d_map                   
        %> 1-d array of bundle types.
        d_mat_map
        %> List of neighbors of each bundle.
        d_neighbors             
        %> Number of bundles.
        d_number_nodes
        %> 1-d array of the number of neigbors each bundle has.
        d_number_neighbors
        %> Bundle width.
        d_h
        %> 1-d array of bundle kinf's.
        d_kinf
        %> K-eigenvalue
        d_k
        %> Fission density
        d_s
        %> 2-d array of fission density (plots, etc.)
        d_f
        %> 1-d array of peaking factors
        d_p
        %> Maximum peaking factor
        d_max_p
        %> Tolerance on eigenvalue.
        d_tolerance_k           
        %> Tolerance on fission source.
        d_tolerance_s           
        % Maximum iterations.
        d_max_iters             
        % Probability of self-collision
        d_wpp                   
        % Probability of leak in a single neighbor
        d_wqp                   
        % Total leakage probability (includes albedo effect)
        d_wleak              
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %> @param   input   User input
        %> @param   mat     Material database
        %> @param   map     Map of assembly type placement
        %> @return          Instance of the FLARE class.
        % ======================================================================
        function this = FLARE(input, mat, map)
            this.d_input = input;
            this.d_h = get(input, 'node_width');
            this.d_mat   = mat;
            number_per_row = zeros(length(map(:, 1)),1);
            number_per_col = number_per_row;
            % Loop through the core map and build the node data
            n = 0; % node counter
            for i = 1:length(map(:, 1))
                for j = 1:length(map(1, :))
                    if (map(i, j) > 0)
                    n = n + 1;
                    % assign my material index
                    this.d_mat_map(n) = map(i, j);
                    % set the map entry to my node index
                    map(i, j) = n;
                    % add one to the number in this row
                    number_per_row(i) = number_per_row(i) + 1;
                    number_per_col(j) = number_per_col(j) + 1;
                    end
                end
            end
            this.d_map = map;
            % subtract first element
            number_per_row(1) = number_per_row(1) - 1;
            % my neighbors
            neighbors = zeros(n, 4);
            % Central node is always adjacent to node 2 all four times.
            neighbors(1, :) = 2;
            % do other nodes
            k = 1;
            for i = 1:length(map(1, :))
                for j = 1:length(map(1, :))
                    if ((i == 1) && (j == 1)) || map(i, j) == 0
                        continue;
                    end
                    k = k + 1;
                    % right
                    if (j < number_per_row(i) + 1)
                        % my right neighbor
                        neighbors(k, 2) = this.d_map(i, j+1);
                    else
                        neighbors(k, 2) = 0; % reflector
                    end
                    % left
                    if j > 2
                        neighbors(k, 1) = this.d_map(i, j-1);
                    else
                        % rotational
                        neighbors(k, 1) = this.d_map(1, i);
                    end
                    % bottom
                    if (i < number_per_col(j) + 1)
                        % my bottom neighbor
                        neighbors(k, 3) = this.d_map(i+1, j);
                    else
                        neighbors(k, 3) = 0; % reflector
                    end
                    % top
                    if (i > 1)
                        neighbors(k, 4) = this.d_map(i-1, j);
                    else
                        neighbors(k, 4) = this.d_map(j, 2);
                    end 
                end
            end
            this.d_neighbors = neighbors;
            this.d_number_nodes = n;
            % Now, compute the coefficients, wpp and wqp.  
            g   = get(input, 'mixing_factor');
            aI  = get(input, 'albedo_single');
            aII = get(input, 'albedo_double');
            this.d_wleak = zeros(n, 1);
            this.d_wpp   = zeros(n, 1);
            this.d_wqp   = zeros(n, 1);
            this.d_kinf  = zeros(n, 1);
            for i = 1:n
                % Compute the migration area
                id   = this.d_mat_map(i);
                D1   = diff_coef(mat, id, 1);
                D2   = diff_coef(mat, id, 2);
                sig1 = sigma_t(mat, id, 1) - sigma_s(mat, id, 1, 1);
                sig2 = sigma_t(mat, id, 2) - sigma_s(mat, id, 2, 2);
                M = sqrt(D1/sig1 + D2/sig2);
                this.d_kinf(i) = (nu_sigma_f(mat, id, 1) + ...
                                  nu_sigma_f(mat, id, 2) * ...
                                  sigma_s(mat, id, 2, 1)  / sig2) / sig1;
                % Compute a single probability for me
                w = (1 - g) * 0.5 * M / this.d_h + g * M^2 / this.d_h^2;
                % Compute the this-probabilities
                this.d_number_neighbors(i) = sum(neighbors(i, :)>0);
                % Subtract leakage to actual neigbors
                this.d_wpp(i)    = 1.0 - this.d_number_neighbors(i) * w; 
                if 4 - this.d_number_neighbors(i) == 1
                    this.d_wpp(i) =this.d_wpp(i) - w*(1-aI);
                    this.d_wleak(i)  = w*(1-aI);
                elseif 4 - this.d_number_neighbors(i) == 2
                    this.d_wpp(i) = this.d_wpp(i) - 2*w*(1-aII);
                    this.d_wleak(i)  = 2*w*(1-aII);
                elseif this.d_number_neighbors(i) == 4
                    this.d_wleak(i)  = 0.0;
                else
                    error('number of neighbors is off!!')
                end
                % My leakage to another
                this.d_wqp(i) = w;  
            end
           
        end
        
        % ======================================================================
        %> @brief Solve the nodal equations
        % ======================================================================
        function [s, k] = solve(this)
           
            % Initialize the fission source and normalize
            s = ones(this.d_number_nodes, 1);
            s = s / norm(s);
            
            % Guess k = 1
            k = 1.0;
            
            % Outer iteration
            for j = 1:100
                
                % Jacobi, sInner iteration (just one)
                s_o = s;
                for p = 1:this.d_number_nodes
                    s(p) = this.d_wpp(p) * s_o(p);
                    for q = 1:4
                        if this.d_neighbors(p, q) == 0
                            continue
                        end
                        qq = this.d_neighbors(p, q);
                        s(p) = s(p) + this.d_wqp(qq)*s_o(qq)* this.d_kinf(p)/k;
                    end
                    s(p) = s(p) * this.d_kinf(p)/k;
                end
                
                % Update k
                k_o = k;
                k   = (sum(s) - sum(s.*this.d_wleak)) / sum(s./this.d_kinf);
                k2  = k_o*norm(s)/norm(s_o);
                s   = s / norm(s);
                kerr = abs(k-k_o);
                serr = norm(s-s_o); 
%                 disp([' it = ', num2str(j), ' k = ', num2str(k), ...
%                       ' kerr = ', num2str(kerr), ' k2 = ', num2str(k2), ...
%                       ' serr = ', num2str(serr)  ])
                if  serr < 1e-4
                    break
                end
            end
            this.d_s = s;
            this.d_k = k;
        end
        
        function make_power(this)
            f = this.d_map * 0.0;
            for i = 1:length(f(:, 1));
                for j = 1:length(f(1, :));
                    if this.d_map(i, j) == 0
                        continue
                    else
                       f(i, j) = this.d_s(this.d_map(i, j)); 
                    end
                end
            end
            this.d_f = f;
            this.d_f(:, 1) = this.d_f(1, :);
        end
        
        function plot_peak(this)
            make_power(this);
            mean_s = mean(this.d_s);
            this.d_p = this.d_s / mean_s;
            this.d_max_p = max(this.d_p);
            pcolor(this.d_f/mean_s)
            axis square
        end
        
    end
    
end
