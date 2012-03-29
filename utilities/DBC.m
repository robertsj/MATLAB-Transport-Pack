%> @file  DBC.m
%> @brief DBC class definition.
% ==============================================================================
%> @brief Provides some basic Design-By-Contract functionality.
%
%> Note, these function calls can be expensive.  I'm not sure of the best
%> way to turn them "off".  A return statement placed before the evaluation
%> helps, but better would be to comment the DBC call.  Perhaps a
%> preprocessor would be good.
% ==============================================================================
classdef DBC
   
    methods (Static)
        
        function Require(testCondition)
            if (~evalin('caller', testCondition))
                error(['REQUIRE FAILED: ' testCondition]);
            end
        end
        
        function Ensure(testCondition)
            if (~evalin('caller', testCondition))
                error(['ENSURE FAILED: ' testCondition]);
            end
        end
        
        function Assert(testCondition)
            if (~evalin('caller', testCondition))
                error(['ASSERT FAILED: ' testCondition]);
            end
        end
        
        function Insist(testCondition, message)
            if (~evalin('caller', testCondition))
                error(['INSIST FAILED: ',testCondition, ' ---> ', message]);
            end
        end
        
    end
    
end
