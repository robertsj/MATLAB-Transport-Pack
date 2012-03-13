%> @file  ERME_Solver.m
%> @brief ERME_Solver class definition.
% ==============================================================================
%> @brief Abstract solver class for eigenvalue response matrix equations.
% ==============================================================================
classdef ERME_Solver < handle
    
    properties (Access = protected)
        %> Input database
        d_input
        %> ERME problem
        d_problem  
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
        %> @return              Instance of the Connect class.
        % ======================================================================
        function this = ERME_Solver(input, problem)
            this.d_input   = input;
            this.d_problem = problem;
            
           
        end
        
        % ======================================================================
        %> @brief  Solve the problem.
        % ======================================================================        
        this = solve(this);
        
        % ======================================================================
        %> @brief  Initialize current (normalized uniform zeroth order)
        % ====================================================================== 
        function J = init_J(this)
            J = zeros(size_J(this.d_problem), 1);
            ng = get(this.d_input, 'number_groups');
            % Degrees of freedom within a group on a surface
            dof_surface = size_J(this.d_problem)             / ...
                          number_elements(this.d_problem)    / ...
                          ng / ...
                          number_faces(this.d_problem); 
            % This assumes that the group is the outer variable
            for g = 1:ng;
                J(g:(dof_surface)*ng:end, 1) = 1.0;
            end 
            J = J / norm(J);
        end
        
        % ======================================================================
        %> @brief  Compute the nonlinear residual.  
        %
        %> For simplicity, this can be used in convergence criteria in Picard
        %> and Newton iterations.  It computes
        %> @f[
        %>     \mathbf{f(x)} =  
        %>       \left [\begin{array}{c}
        %>         (\mathbf{M}\mathbf{R}(k)-\lambda \mathbf{I}) \mathbf{J_-} \\
        %>         \mathbf{F}(k)\mathbf{J_-} - (k\mathbf{L}(k)\mathbf{J_-} ) \\
        %>         \frac{1}{2} \mathbf{J^T_-} \mathbf{J_-} - \frac{1}{2}  
        %>       \end{array} \right ]  = 
        %>     \mathbf{0} \, ,
        %> @f]
        %>
        % ====================================================================== 
        function f = residual(this, x)
            % Decompose unknown vector for clarity.
            J = x(1:end-2, 1);
            k = x(end-1);
            lambda = x(end); 
            % Update the operators (won't if k is the same!) and get them.
            update(this.d_problem, k);
            [R, F, A, L, M, leak] = get_operators(this.d_problem);
            % Compute the residual 
            f_J = M*(R*J) - lambda*J;
            f_k = F*J - k*(A*J + leak*(L*J));
            f_lambda = 0.5 - 0.5*(J'*J);
            f = [f_J; f_k; f_lambda];
        end
        
    end
    
end