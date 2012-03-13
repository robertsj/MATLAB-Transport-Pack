%> @file  ERME_Picard.m
%> @brief ERME_Picard class definition.
% ==============================================================================
%> @brief Solve eigenvalue response matrix equations via Picard iteration.
%
%> Currently, only acceleration via Steffensen extrapolation is implemented.
%> Other ideas being investigated are k updated via a Rayleigh quotient, 
%> CMFD/p-CMFD, and low order respons
% ==============================================================================
classdef ERME_Picard < ERME_Solver
    
    properties
        %> Residual history
        d_norm_residual_hist
        %> Final solution
        d_J
        d_k
        d_lambda
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief  Class constructor
        %
        %> Is is *assumed* that the element map is corrected.  The .m file based
        %> input has defined maps in real orientation; this has to be corrected
        %> during parsing.
        %>
        %> @param  input        Input database.
        %> @param  elements     Element map.
        %> @return              Instance of the ERME_Picard class.
        % ======================================================================
        function this = ERME_Picard(input, problem)
            this = this@ERME_Solver(input, problem);
        end
        
        % ======================================================================
        %> @brief  Solve the problem.
        % ======================================================================
        function this = solve(this)
            
            % Set convergence criteria.
            inner_tolerance = get(this.d_input, 'inner_tolerance');
            inner_max_iters = get(this.d_input, 'inner_max_iters');
            outer_tolerance = get(this.d_input, 'outer_tolerance');
            outer_max_iters = get(this.d_input, 'outer_max_iters'); 
            
            % Initialize current.
            J = init_J(this);
            
            % Initialize k-eigenvalue and lambda-eigenvalue
            k = 1.0;
            lambda = 1.0;
            
            % Get response matrices for initial k
            update(this.d_problem, k);
            [R, F, A, L, M, leak] = get_operators(this.d_problem);
            
            % Get initial nonlinear residual norm.
            norm_residual = norm(residual(this, [J; k; lambda]));
            % and start a residual history.
            norm_residual_hist = zeros(outer_max_iters + 1, 1);
            norm_residual_hist(1) = norm_residual;
            
           
            
            % ==================================================================
            % Outer iterations
            % ==================================================================
            
            iteration = 0;
            
            while iteration < outer_max_iters && norm_residual > outer_tolerance
                
                fprintf(1,'it = %d, keff = %12.10f, lambda = %12.10f, norm = %12.10e\n', ...
                    iteration, k, lambda, norm_residual);
                
                iteration = iteration + 1;
                
                % ==============================================================
                % Inner iterations
                % ==============================================================
                MR = M*R; % Set a single operator.                         
                opts.disp  = 0;
                opts.tol   = inner_tolerance;
                opts.maxit = inner_max_iters;
                [J, lambda] = eigs(MR, 1, 'LM', opts);                    
                J = J*sign(sum(J)); % We want the positive direction                           

                % ==============================================================
                % Eigenvalue update
                % ==============================================================
                %k_old       = k;                              % store keff
                %J_old       = J;                              % store currents
                %lambda_old  = lambda;                         % store current eigenvalue
                gain        = F*J;                            % compute gains
                absorb      = A*J;                            % absorption
                leak        = leak*(L*J);                     % leakage
                loss        = leak + absorb;                  % total loss
                k           = gain / loss;                    

                % Get response matrices for initial k
                update(this.d_problem, k);
                [R, F, A, L, M, leak] = get_operators(this.d_problem);
                
                % Get initial nonlinear residual norm.
                norm_residual = norm(residual(this, [J; k; lambda]));
                % and start a residual history.
                norm_residual_hist(iteration+1) = norm_residual;
                
            end
            
            disp('*** final results ***')
            fprintf(1,'it = %d, keff = %12.10f, lambda = %12.10f, norm = %12.10e\n', ...
                iteration, k, lambda, norm_residual);
            
            % Store only those needed.
            this.d_norm_residual_hist = norm_residual_hist(1:iteration+1);
            
            % Store the results
            this.d_J = J;
            this.d_k = k;
            this.d_lambda = lambda;

        end
        
    end
    
end
% 
% template <class Inner>
% scalar PowerIter<Inner>::Aitken(scalar k0, scalar k1, scalar k2)
% {
%   scalar kA = k0 - pow(k1 - k0, 2) / (k2 - 2.0 * k1 + k0);
%   return kA;
% }

