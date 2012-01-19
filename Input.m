%> @file  Input.m
%> @brief Input class definition.
% ==============================================================================
%> @brief Stores and processes user input from a driver script.
%
%> Finish me.
% ==============================================================================
classdef Input < handle
    
    properties (Access = private)
        %> Input data map, consisting of "keys" with "values"
        d_map = containers.Map();
    end
         
    % Input parameters.  Some of these are optional/defaulted, as denoted.
    properties (Access = public)
        
        % Boundary conditions.  The available options are vacuum and reflected
        % for general problems.  For response function generation ...
        % 
        bc_left   = 'vacuum';
        bc_right  = 'vacuum';
        bc_bottom = 'vacuum';
        bc_top    = 'vacuum';
        
        % Materials
        number_groups    = 1;
        
        % Discretization
        equations = 'DD'; % DD, SD, SC, or varied
        
        % Convergence criteria.
        inner_max_iters  = 900;     % Maximum inner (within-group) iterations
        inner_tolerance  = 1e-5;    % Pointwise flux relative error tolerance
        max_outer        = 100;     % Maximum outer (group) iterations. 
        tol_outer        = 1e-5;    % Pointwise flux relative error tolerance
        max_fission      = 100;     % Maximum fission source iterations.  This 
                                    %   applies to eigenproblems and fixed
                                    %   source problems with multiplication.
        tol_fission      = 1e-3;    % Pointwise fission source relative error 
                                    %   tolerance
        tol_eigenvalue   = 1e-4;    % Absolute tolerance on eigenvalue.
        
        % Quadrature selection
        quadrature_type = QuadratureTypes.LEVELSYMMETRIC; % Quadrature type.
        quadrature_order = 8;                             % Quadrature order.
        
        % Problem type
        problem_type = 'Fixed'           % Fixed, FixedMult, Eigen, Response
        solver_type  = 'SourceIteration'; % Only option for now.
        
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> @return Instance of the Mesh class with standard defaults.
        % ======================================================================
        function obj = Input()
            
            % Materials
            put(obj, 'number_groups',   1);
            put(obj, 'downscatter',     0);
            
            % Boundary conditions.
            put(obj, 'bc_left',         'vacuum');
            put(obj, 'bc_right',        'vacuum');
            put(obj, 'bc_top',          'vacuum');
            put(obj, 'bc_bottom',       'vacuum');
            
            % Discretization
            put(obj, 'equations',       'DD');
            put(obj, 'dimension',       2);
            
            % Quadrature
            put(obj, 'quadrature_type',     'LevelSymmetric');
            put(obj, 'quadrature_order',    8);
            
            % Solver 
            put(obj, 'problem_type',        'Fixed');  
            
            % Inner solver
            put(obj, 'inner_solver',        'SI'); 
            put(obj, 'inner_tolerance',     1e-6); 
            put(obj, 'inner_max_iters',     1000); 
        end
        
        % ======================================================================
        %> @brief Put a new (key, value) pair into the input database.
        %
        %> @return Instance of the Mesh class with standard defaults.
        % ======================================================================
        function obj = put(obj, key, value)
        	obj.d_map(key) = value;
        end
        
        % ======================================================================
        %> @brief Get the value for a key in the input database
        %
        %> @return Value for the key or zero along with a warning.
        % ======================================================================
        function value = get(obj, key)
            if (isKey(obj.d_map, key))
                value = obj.d_map(key); 
            else
                value = 0;
            end
        end
        
        % ======================================================================
        %> @brief Check if a key exists.
        %
        %> @return  1 if it exists, 0 otherwise.
        % ======================================================================
        function value = contains(obj, key)
            if (isKey(obj.d_map, key))
                value = 1; 
            else
                value = 0;
            end
        end
        
    end
    
end