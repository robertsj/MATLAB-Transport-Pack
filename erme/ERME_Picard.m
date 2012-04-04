%> @file  ERME_Picard.m
%> @brief ERME_Picard class definition.
% ==============================================================================
%> @brief Solve eigenvalue response matrix equations via Picard iteration.
%
%> Currently, only acceleration via Steffensen extrapolation is implemented.
%> Other ideas being investigated are k updated via a Rayleigh quotient, 
%> CMFD/p-CMFD, and low order response
%>
%> Input database options of relevance include
%>   - erme_steffensen   (1 turns it on)
%>   - erme_picard_inner (One of 'eigs', 'krylovschur', 'arnoldi', 'power')
%> plus possibly others general to @ref ERME_Solver.
%>
% ==============================================================================
classdef ERME_Picard < ERME_Solver
    
    properties
        %> Residual history
        d_norm_residual_hist
        %> 
        d_number_iterations
        %> Inner solver (eigs, krylovschur, arnoldi, power)
        d_inner_solver
        %> SLEPc EPS object
        d_eps
        %> PETSc shell matrix for MR
        d_MR
        %> PETSc Vec for J
        d_pJ
        %> Do we use Steffensen acceleration?
        d_steffensen
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
        %> @param  problem      ERME object.
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
            this.d_steffensen      = get(this.d_input, 'erme_steffensen'); 
            this.d_inner_solver    = get(this.d_input, 'erme_picard_inner');
            
            if ~this.d_inner_solver
                this.d_inner_solver = 'eigs';
            end    

            % Initialize current.
            J = init_J(this);
            
            % Initialize k-eigenvalue and lambda-eigenvalue
            k =  get(this.d_input, 'erme_initial_keff'); 
            if ~k
                k = 1.0;
            end
            k_1 = 0; % keff one iteration ago
            k_2 = 0; % keff two iterations ago
            lambda = 1.0;
            
            % Get response matrices for initial k
            update(this.d_problem, k);
            [R, F, A, L, M, leak] = get_operators(this.d_problem);
            
            % Get initial nonlinear residual norm.
            norm_residual = norm(residual(this, [J; k; lambda]));
            % and start a residual history.
            norm_residual_hist = zeros(this.d_outer_max_iters + 1, 1);
            norm_residual_hist(1) = norm_residual;
            
            % If using SLEPc, initialize
            use_slepc = 0;
            if ~strcmp(this.d_inner_solver, 'eigs')
                use_slepc = 1;
                % Setup SLEPc solver.
                this.setup_slepc();
            end
            
            % ==================================================================
            % Outer iterations
            % ==================================================================
            
            iteration = 0;

            while iteration < this.d_outer_max_iters && ...
                  norm_residual > this.d_outer_tolerance
                
                fprintf(1, ['it = %d, keff = %12.10f, lambda = %12.10f,', ...
                    ' norm = %12.10e\n'], iteration, k, lambda, norm_residual);
                
                iteration = iteration + 1;
                
                % ==============================================================
                % Inner iterations
                % ==============================================================
       
                
                if ~use_slepc
                    % Using eigs
                    opts.disp  = 0;
                    opts.tol   = this.d_inner_tolerance;
                    opts.maxit = this.d_inner_max_iters;
                    MR = M*R; % Set a single operator.  
                    [J, lambda] = eigs(MR, 1, 'LM', opts);     
                    
                else
                    % Using SLEPc
                    % Remake operator.
                    this.d_MR(:, :) = M*R; 
                    % Solve and extract solution.
                    this.d_eps.Solve();
                    lambda = this.d_eps.GetEigenpair(1, this.d_pJ);
                    J = this.d_pJ(:);
                    % Destroy operator and solver.
                    %this.d_eps.Destroy();
                    %this.d_MR.Destroy();
                end
                J = J*sign(sum(J)); % We want the positive direction  
                
                % ==============================================================
                % Eigenvalue update
                % ==============================================================
                k_2         = k_1;                            % Keep the last
                k_1         = k;                              %   two keffs
                gain        = F*J;                            % compute gains
                absorb      = A*J;                            % absorption
                leak        = leak*(L*J);                     % leakage
                loss        = leak + absorb;                  % total loss
                k           = gain / loss;                    

                % Steffensen update if required.
                if this.d_steffensen && iteration > 2
                   k = this.extrapolate(k, k_1, k_2);
                end

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

            this.d_problem.update_state(J, k, lambda);
            this.d_number_iterations = iteration;
            
            % Destroy SLEPc stuff if applicable
            if use_slepc
                this.d_pJ.Destroy();
                this.d_MR.Destroy();
                this.d_eps.Destroy();
                SlepcFinalize();
            end
            
        end
        
    end
    
    methods (Access = private)
        
        % ======================================================================
        %> @brief  Steffensen/Aitken extrapolation
        %> @param   k       Computed eigenvalue
        %> @param   k_1     Eigenvalue computed 1 time ago
        %> @param   k_2     Eigenvalue computed 2 times ago
        %> @return  k_ex    Extrapolated eigenvalue
        % ======================================================================   
        function k_ex = extrapolate(this, k, k_1, k_2)  
            k_ex = k_2 - (k_1 - k_2)^2 / (k - 2.0 * k_1 + k_2);
        end
        
        
        % ======================================================================
        %> @brief  Initialize SLEPc objects for inner eigensolve.
        % ======================================================================   
        function setup_slepc(this)        
            disp('Setting up SLEPc')
            % Initialize SLEPc
            slepc_options = get(this.d_input, 'slepc_options');
            SlepcInitialize(slepc_options);
            % Initialize operator
            [R, F, A, L, M, leak] = get_operators(this.d_problem);
            this.d_MR = PetscMat(M*R);
            % Make new solver.
            this.d_eps = SlepcEPS();
            % Non-hermitian.
            this.d_eps.SetProblemType(SlepcEPS.NHEP);
            % Number
            this.d_eps.SetDimensions(1, 3);
            % Largest magnitute.
            this.d_eps.SetWhichEigenpairs(1);
            % Set type.
            this.d_eps.SetType(this.d_inner_solver);
            % Set convergence criteria, operator, and any options.
            this.d_eps.SetTolerances(this.d_inner_tolerance, ...
                                     this.d_inner_max_iters);
            this.d_eps.SetOperators(this.d_MR);
            this.d_eps.SetFromOptions();
            % Create vector for extracting eigenvector.
            n = size_J(this.d_problem);
            this.d_pJ = PetscVec(zeros(n, 1));
        end
        
    end
    
end


