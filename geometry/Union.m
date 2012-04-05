%> @file  Union.m
%> @brief Union class definition.
% ==============================================================================
%> @brief 
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef Union
    
    methods (Static)
        
        function b = evaluate(left, right, location)
            b = left.contains(location) || right.contains(location);
        end
        
    end
    
end
        