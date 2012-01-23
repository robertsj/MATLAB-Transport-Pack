%> @file  SCMOC.m
%> @brief SCMOC class definition.
% ==============================================================================
%> @brief Step characteristic approximation for MOC.
%
%> In the method of characteristics, the flux is solved for along a
%> track assuming a flat source.  For a given incident flux into 
%> a track segment, we can define the outgoing segment flux
%> \f[
%>     \psi_{out} = A\psi_{in}  + B Q   \, ,
%> \f] 
%> and average segment flux
%> \f[
%>     \bar{\psi} = \frac{1}{l} \Big ( B \psi_{in} +  C Q \Big )  \, ,
%> \f] 
%> where 
%> \f[
%>     A = e^{-\Sigma_t \tau} \, ,
%> \f]
%> \f[
%>     B = \frac{1}{\Sigma_t} ( 1- A ) \, ,
%> \f]
%> and
%> \f[
%>     C = \frac{l}{\Sigma_t} \Big( 1- \frac{1-A}{\tau} \Big ) \, ,
%> \f]
%> where \f$ l \f$ is the segment length and 
%> \f$ \tau = \Sigma_t l  \f$ is optical path length.
%>
%> The step characteristic method is positive but only first-order
%> accurate in space.
%>
%> Reference:  A. Hebert, <em>Applied Reactor Physics</em>.
%>
%> \sa DDMOC
%> 
% ==============================================================================
classdef SCMOC < Equation
    
    properties

    end
    
    properties (Constant)

    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        %
        %> Set the mesh and material.
        %>
        %> @param mesh        	Problem mesh.
        %> @param mat          	Material definitions.
        %>
        %> @return Instance of the SCMOC class.
        % ======================================================================
        function obj = SCMOC(mesh, mat)
            % Call the base class.
            obj = obj@Equation(mesh, mat); 
        end
        
        % ======================================================================
        %> @brief Setup the equations for a group.
        %
        %> Here, we'll go through the grid and produce a fine mesh matrix
        %> of total cross-sections.  This isn't the best thing for memory,
        %> but it cuts down a lot on the time within solve.
        %>
        %> @param group     Current group.
        % ======================================================================
        function obj = setup_group(obj, group)
            % Nothing here for now.
        end
        
        % ======================================================================
        %> @brief Setup the equations for an octant.
        %>
        %> @param octant    Current octant.
        % ======================================================================
        function obj = setup_octant(obj, octant)
            % Nothing here for now.   
        end
        
        % ======================================================================
        %> @brief Setup the equations for an angle.
        %>
        %> @param phi   Azimuth with respect to x axis.
        % ======================================================================
        function obj = setup_angle(obj, phi, ~)
            % Get the widths from mesh.
        end
        
        % ======================================================================
        %> @brief Solve for the cell-center and outgoing edge fluxes.
        %>
        %> @param psi_in    Incident flux vector
        %> @param s         Region isotropic source
        %> @param sig       Region total cross-section
        %> @param t         Segment length (includes polar scaling)
        %> @return          Segment exit and average angular flux.
        % ======================================================================
        function [psi_out, psi_avg] = solve(obj, psi_in, s, sig, t)   
            % Many of these values (A,B,C,etc) should be
            % able to be precomputed.
            A    = exp(-sig*t);                 % A, B, and C
            B    = (1.0-A)/sig;                 % from Hebert's book
            C    = t/sig * (1 - (1-A)/(t*sig));
            psi_out = A*psi_in + B*s;           % Segment exiting flux.
            if (psi_out < 0)
                error('negative flux cannot happen with SC!')
            end
            psi_avg = (1/t) * (B*psi_in + C*s);	% Average segment flux.
        end
        
    end
    
end