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
        %> Inner tolerance
        d_inner_tolerance
        %> Maximum inner iterations
        d_inner_max_iters
        %> Outer tolerance
        d_outer_tolerance
        %> Maximum outer iterations
        d_outer_max_iters
        
        %> Final expansion coefficients
        d_J
        %> Final eigenvalue
        d_k
        %> Final current eigenvalue
        d_lambda
        
        d_time = 0;
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
            % Set convergence criteria
            this.d_inner_tolerance = get(this.d_input, 'inner_tolerance');
            if ~this.d_inner_tolerance
                this.d_inner_tolerance = 1e-8;
            end
            this.d_inner_max_iters = get(this.d_input, 'inner_max_iters');
            if ~this.d_inner_max_iters
                this.d_inner_max_iters = 100;
            end
            this.d_outer_tolerance = get(this.d_input, 'outer_tolerance');
            if ~this.d_outer_tolerance
                this.d_outer_tolerance = 1e-8;
            end
            this.d_outer_max_iters = get(this.d_input, 'outer_max_iters'); 
            if ~this.d_outer_max_iters
                this.d_outer_max_iters = 20;
            end
        end
        
        % ======================================================================
        %> @brief  Solve the problem.
        % ======================================================================        
        this = solve(this);
        
        function k = get_keff(this)
            k = this.d_k;
        end
        
        function l = get_lambda(this)
            l = this.d_lambda; 
        end
        
        % ======================================================================
        %> @brief  Initialize unknowns
        % ====================================================================== 
        function [J, k, lambda] = init(this)
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
            % Initialize k-eigenvalue and lambda-eigenvalue
            k =  get(this.d_input, 'erme_initial_keff'); 
            if ~k
                k = 1.0;
            end
            lambda = 1.0;
        end

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

        function p = problem(this)
            p = this.d_problem; 
        end
        
        function t = mytime(this)
            t = this.d_time; 
        end
        
    end
    
end
