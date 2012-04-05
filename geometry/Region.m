%> @file  Region.m
%> @brief Region class definition.
% ==============================================================================
%> @brief Represents a single flat source region.
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef Region < handle
    
    properties (Access = public)
        %> name of the region
        name
        %> actual volume of the region
        volume = 1;
        %> calculated volume based on segments within region
        volume_calc
        %> scalar flux in region
        flux  
        %> source that each track picks up in this region
        source 
        %> cell array of segments within region
        segments 
        %> @ref SurfaceNode @ref OperatorNode representing a closed volume
        node
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the Geometry class.
        % ======================================================================
        function self = Region(name)
            self.name           = name;
            self.volume         = 1.;
            self.volume_calc    = 0.;
            self.flux           = 0.;
            self.source         = 0.;
            self.segments       = {};
            self.node           = [];
        end
        
    end % public methods
    
end
            
            