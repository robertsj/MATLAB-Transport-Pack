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
            path(path,'/home/robertsj/opt/petsc/petsc-3.2-p5/bin/matlab/classes/')
            PetscInitialize({'-snes_monitor','-snes_ls', 'basic','-pc_type','none'});
            
            inner_tolerance = get(this.d_input, 'inner_tolerance');
            inner_max_iters = get(this.d_input, 'inner_max_iters');
            outer_tolerance = get(this.d_input, 'outer_tolerance');
            outer_max_iters = get(this.d_input, 'outer_max_iters'); 

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
            
            % Seed with one inner.
            MR = M*R;                     
            opts.disp  = 0;
            opts.tol   = inner_tolerance;
            opts.maxit = inner_max_iters;
            [J, lambda] = eigs(MR, 1, 'LM', opts);                    
            J = J*sign(sum(J)); % We want the positive direction  
            gain        = F*J;                            % compute gains
            absorb      = A*J;                            % absorption
            leak        = leak*(L*J);                     % leakage
            loss        = leak + absorb;                  % total loss
            k           = gain / loss;                    

            % Set the unknown
            x = [J;k;lambda];

            % ==================================================================
            % PETSc SNES
            % ==================================================================

            n = length(x);

            % Create work vector for nonlinear solver and location for solution
            v = PetscVec();
            v.SetType('seq');
            v.SetSizes(n);
            z = v.Duplicate();
            disp('lala')
            z = PetscVec(x);
           % z.SetValues(1:n);
           % sum(z(1:10))

            % Create a matrix for the Jacobian for Newton method
            Jac = PetscMat();
            Jac.SetType('seqaij');
            Jac.SetSizes(n, n);
            Pre = PetscMat(speye(n));
%            Pre.SetType('seqaij');
%            Pre.SetSizes(n, n);
            
            % Create the nonlinear solver (line search)
            snes = PetscSNES();
            snes.SetType('ls');
            
            %  Provide a function 
            snes.SetFunction(v, 'nonlinear_residual_petsc', this);

            %  Provide a function that evaluates the Jacobian
            snes.SetJacobian(Jac, Jac, 'jacobian_petsc', this);

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
            v.Destroy();
            z.Destroy();
            Jac.Destroy();
            Pre.Destroy();
            snes.Destroy();
            PetscFinalize();
        end

        
    end
    
end


