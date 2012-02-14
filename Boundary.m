%> @file  Boundary.m
%> @brief Boundary class definition.
% ==============================================================================
%> @brief Holds the incident and outgoing boundary fluxes.
%
%> More here...
% ==============================================================================
classdef Boundary < handle
    
    properties (Constant)
       IN   =  1;
       OUT  = -1;
    end

    properties (Access = protected)
        %> Input
        d_input
        %> Spatial mesh.
        d_mesh
        %> Angular mesh.
        d_quadrature
        %> Boundary flux cell array
        d_boundary_flux
        %> Cell array of \ref BoundaryCondition thisects for each surface
        d_bc  
        %> Current group.
        d_g = 0;
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %>
        %> All boundary flux arrays are sized, and the boundary conditions are
        %> constructed for each surface.
        %>
        %> @param input         User input.
        %> @param mesh          Spatial mesh.
        %> @param quadrature  	Angular mesh.
        %>
        %> @return Instance of the Boundary class.
        % ======================================================================
        function this = Boundary(input, mesh, quadrature)
        
            this.d_input         = input;
            this.d_mesh          = mesh;
            this.d_quadrature    = quadrature;

            % Boundary conditions
            this.d_bc = cell(2*mesh.DIM);
            
            if strcmp(get(input, 'bc_left'), 'vacuum')
                this.d_bc{mesh.LEFT}   = ...
                    Vacuum(this, input, mesh, quadrature, mesh.LEFT);
            elseif strcmp(get(input, 'bc_left'), 'reflect')
                this.d_bc{mesh.LEFT}   = ...
                    Reflective(this, input, mesh, quadrature, mesh.LEFT);
            elseif strcmp(get(input, 'bc_left'), 'response')
                this.d_bc{mesh.LEFT}   = ...
                    Response(this, input, mesh, quadrature, mesh.LEFT);
            end
            if strcmp(get(input, 'bc_right'), 'vacuum')
                this.d_bc{mesh.RIGHT}   = ...
                    Vacuum(this, input,  mesh, quadrature, mesh.RIGHT);
            elseif strcmp(get(input, 'bc_right'), 'reflect')
                this.d_bc{mesh.RIGHT}   = ...
                    Reflective(this, input, mesh, quadrature, mesh.RIGHT);
            end   
            if mesh.DIM > 1
                if strcmp(get(input, 'bc_bottom'), 'vacuum')
                    this.d_bc{mesh.BOTTOM}   = ...
                        Vacuum(this, input, mesh, quadrature, mesh.BOTTOM);
                elseif strcmp(get(input, 'bc_bottom'), 'reflect')
                    this.d_bc{mesh.BOTTOM}   = ...
                        Reflective(this, input, mesh, quadrature, mesh.BOTTOM);
                end
                if strcmp(get(input, 'bc_top'), 'vacuum')
                    this.d_bc{mesh.TOP}   = ...
                        Vacuum(this, input, mesh, quadrature, mesh.TOP);
                elseif strcmp(get(input, 'bc_top'), 'reflect')
                    this.d_bc{mesh.TOP}   = ...
                        Reflective(this, input, mesh, quadrature, mesh.TOP);
                end
            end
            
        end
        
        % ======================================================================
        %> @brief Initialize for a within-group solve.
        %>
        %> More detailed description of what the constructor does.
        %>
        %> @param input         User input.
        %> @param mesh          Spatial mesh.
        %> @param quadrature  	Angular mesh.
        %>
        %> @return Instance of the Boundary class.
        % ======================================================================
        function this = initialize(this, g)
            this.d_g = g;
        end
        
        
        % ======================================================================
        %> @brief Set the incident boundary fluxes.
        %>
        %> This is called after every sweep.  The sweeper updates all the
        %> outgoing angular fluxes.  This routine updates the incident fluxes
        %> according to the boundary condition.
        %>
        %> @param order         Quadrature order. This differs from quadrature 
        %>                      to quadrature, but e.g. for level symmetric, 
        %>                      it's the number of unique directional cosines.
        %>
        %> @return Instance of the Quadrature class.
        % ======================================================================  
        function this = set(this)
            for side = 1:2*this.d_mesh.DIM;
                update(this.d_bc{side});
            end
        end
        
        % ======================================================================
        %> @name Getters
        %> @{
        % ======================================================================
        
        % ======================================================================
        %> @brief Get a vertical face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch 
        %> @return      Boundary flux array.
        % ======================================================================
        f = get_psi_v(this, o, a, inout)
        
        % ======================================================================
        %> @brief Get a horizontal face incident/exiting boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param inout Incoming/outgoing switch         
        %> @return      Boundary flux array.
        % ======================================================================
        f = get_psi_h(this, o, a, inout)
        
        % @}
        
        % ======================================================================
        %> @name Setters
        %> @{
        % ======================================================================
        
        % ======================================================================
        %> @brief Set a vertical face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array.
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        % example: sweeping left to right, i.e octant 1.  At the end, we set the
        % right boundary octant 1 flux.  inout=OUT=-1.  octant(1) = 1.  Hence,
        % the logic below finds RIGHT.  
        %
        this = set_psi_v(this, o, a, f, inout)
        
        % ======================================================================
        %> @brief Set a horizontal face incident boundary flux.
        %
        %> @param o     Octant index.
        %> @param a     Angle index (within an octant).
        %> @param f     Flux array.
        %> @param inout Incoming/outgoing switch         
        % ======================================================================
        this = set_psi_h(this, o, a, f, inout) 
        
        % @}
    end
 
end