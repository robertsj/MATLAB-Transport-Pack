%> @file  DDMOC.m
%> @brief DDMOC class definition.
% ==============================================================================
%> @brief Diamond difference approximation for MOC.
%
%> In the method of characteristics, the flux is solved for along a
%> track assuming a flat source.  In the diamond difference approximation,
%> the average segment flux is taken as an average of the incident and
%> and exiting flux.  Similar to the step characteristic (SCMOD), we define
%> the outgoing segment flux
%> \f[
%>     \psi_{out} = A\psi_{in}  + B Q   \, ,
%> \f] 
%> and average segment flux
%> \f[
%>     \bar{\psi} = \frac{1}{l} \Big ( B \psi_{in} +  C Q \Big )  \, ,
%> \f] 
%> but now
%> \f[
%>     A = \frac{2-\tau}{2+\tau} \, ,
%> \f]
%> \f[
%>     B = \frac{2l}{2+\tau} \, ,
%> \f]
%> and
%> \f[
%>     C = \frac{l^2}{2+\tau} \, ,
%> \f]
%> where \f$ l \f$ is the segment length and 
%> \f$ \tau = \Sigma_t l  \f$ is optical path length.
%>
%> Unlike the step characteristic approximation, the diamond difference
%> approximation is second order in space but is not guaranteed
%> to yield positive fluxes.
%>
%> Reference:  A. Hebert, <em>Applied Reactor Physics</em>.
%>
%> \sa DDMOC
%> 
% ==============================================================================
classdef DDMOC < Equation
    
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
        %> @return Instance of the DDMOC class.
        % ======================================================================
        function obj = DDMOC(mesh, mat)
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
        function [psi_out, psi_avg] = solve(obj, psi_in, s, sig, t, ~)   
            % Many of these values (A,B,C,etc) should be
            % able to be precomputed.
            tau = sig*t;
            two_p_tau = 2+tau;
            two_m_tau = 2-tau;
            A = two_m_tau / two_p_tau;
            B = 2*l / two_p_tau;            
            C = l^2 / two_p_tau;
            psi_out = A*psi_in + B*s;           % Segment exiting flux.
            if (psi_out < 0)
                warning('negative flux in DD!')
            end
            psi_avg = 0.5*(psi_in + psi_out);	% Average segment flux.
        end
        
    end
    
end