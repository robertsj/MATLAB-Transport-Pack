%> @file  Intersection.m
%> @brief Intersection class definition.
% ==============================================================================
%> @brief 
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef Intersection
    
    methods (Static)
        
        function b = evaluate(left, right, location)
            b = left.contains(location) .* right.contains(location);
        end
        
    end
    
end
        