%> @file  OperatorNode.m
%> @brief OperatorNode class definition.
% ==============================================================================
%> @brief 
%
%> Nearly verbatim translation to MATLAB from P. Romano's Python code.
% ==============================================================================
classdef OperatorNode < handle

    properties
       leftNode
       rightNode
       operator
    end
    
    methods (Access = public)
        
        % ======================================================================
        %> @brief Class constructor
        %> @return Instance of the OperatorNode class.
        % ======================================================================
        function self = OperatorNode(leftNode, rightNode, operator)
            self.leftNode = leftNode;
            self.rightNode = rightNode;
            self.operator = operator;
        end
        
        function b = contains(self, location)
            b = self.operator.evaluate(self.leftNode, self.rightNode, location);
        end
        
    end % public methods
    
end
        