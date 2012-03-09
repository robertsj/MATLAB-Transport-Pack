%> @file  ResponseDriver.m
%> @brief Drives generation of response functions.
% ==============================================================================
%> @brief ResponseDriver class definition.
%
%> This drives generation of response functions for a given set of local
%> problems (aka "elements" or "nodes").  The user specifies a set of values
%> for keff and the maximum orders in space and angle.  The output is placed
%> into an HDF5 file.
%> 
%> 
% ==============================================================================
classdef ResponseDriver < handle

    properties (Access = private)
        d_input
        d_mat
        d_mesh_array
    end
    
    methods
        
        % ======================================================================
        %> @brief Class constructor
        % ======================================================================
        function this = ResponseDriver(input, mat, mesh_array)
            this.d_input = input;
            this.d_mat = mat;
            this.d_mesh_array = mesh_array;
            
            
        end
        
        % ======================================================================
        %> @brief Run
        % ======================================================================
        function this = run(this)

        end
        
    end
    
end