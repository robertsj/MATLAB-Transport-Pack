%> @file  ERME_Newton.m
%> @brief ERME_Newton class definition.
% ==============================================================================
%> @brief Solve eigenvalue response matrix equations via Newton iteration.
%
%> Currently, only acceleration via Steffensen extrapolation is implemented.
%> Other ideas being investigated are k updated via a Rayleigh quotient, 
%> CMFD/p-CMFD, and low order response
%>
%> @todo Add Steffensen extrapolation and possible PETSc/SLEPc extension.
%>
% ==============================================================================
classdef ERME_Newton < ERME_Solver
    
    properties
        %> Residual history
        d_norm_residual_hist
        %> Final expansion coefficients
        d_J
        %> Final eigenvalue
        d_k
        %> Final current eigenvalue
        d_lambda
        d_built = 0
        %> Solution PETSc vector.
        d_x
        %> Type of approximate Jacobian for PC.
        d_approximate_jacobian
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
        %> @return              Instance of the ERME_Newton class.
        % ======================================================================
        function this = ERME_Newton(input, problem)
            this = this@ERME_Solver(input, problem);
        end
        
        % ======================================================================
        %> @brief  Solve the problem.
        % ======================================================================
        function this = solve(this)

            petsc_options = get(this.d_input, 'petsc_options');
            PetscInitialize(petsc_options);

            % Initialize unknown.
            [J, k, lambda] = init(this);
            
            % Get response matrices for initial k
            update(this.d_problem, k);
            [R, F, A, L, M, leak] = get_operators(this.d_problem);
            
            % Get initial nonlinear residual norm.
            norm_residual = norm(residual(this, [J; k; lambda]));
            % and start a residual history.
            norm_residual_hist = zeros(outer_max_iters + 1, 1);
            norm_residual_hist(1) = norm_residual;
            
            % Seed with one crude inner.
            MR = M*R;                     
            opts.disp  = 0;
            opts.tol   = 1e-5;
            opts.maxit = inner_max_iters;
            [J, lambda] = eigs(MR, 1, 'LM', opts);                    
            J = J*sign(sum(J)); % We want the positive direction  
            gain        = F*J;                            % compute gains
            absorb      = A*J;                            % absorption
            leak        = leak*(L*J);                     % leakage
            loss        = leak + absorb;                  % total loss
            k           = gain / loss;                    

            % Set the unknown
            x = [J; k; lambda];

            % ==================================================================
            % PETSc SNES
            % ==================================================================

            n = length(x);

            % Create work vector for nonlinear solver and location for solution
            w = PetscVec();
            w.SetType('seq');
            w.SetSizes(n);
            this.d_x = PetscVec(x);

            % Create a matrix for the Jacobian for Newton method
            Jac = PetscMat();
            Jac.SetSizes(n, n);
            Jac.SetType('shell');
            % Setup the shell, setting 'this' as the PETSc context.
            Jac.ShellSetup(this, 'jacobian_matvec');

            % Set the preconditioner.
            
            
            
%            Pre = PetscMat(speye(n));
%            Pre.SetType('seqaij');
%            Pre.SetSizes(n, n);
            
            % Create the nonlinear solver (line search)
            snes = PetscSNES();
            snes.SetType('ls');
            
            %  Provide a function 
            snes.SetFunction(v, 'nonlinear_residual', this);

            %  Provide a function that evaluates the Jacobian
            snes.SetJacobian(Jac, Jac, 'jacobian_update', this);

            %  Solve the nonlinear system
            snes.SetFromOptions();
            snes.Solve(z);
            %x.View;
            x      = z(:);
            J      = x(1:end-2, 1);
            k      = x(end-1);
            lambda = x(end); 
            disp('*** final results ***')
            fprintf(1,'it = %d, keff = %12.10f, lambda = %12.10f, norm = %12.10e\n', ...
                1, k, lambda, norm_residual);
            snes.View;

            % Cleanup PETSc
            v.Destroy();
            z.Destroy();
            Jac.Destroy();
            Pre.Destroy();
            snes.Destroy();
            PetscFinalize();
        end

        
    end
    
end


