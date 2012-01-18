classdef DBC
    % Provides several Design-By-Contract functions
   
    methods (Static)
        
        function Require(testCondition)
            return
            if (~evalin('caller', testCondition))
                error(['REQUIRE FAILED: ' testCondition]);
            end
        end
        
        function Ensure(testCondition)
            return
            if (~evalin('caller', testCondition))
                error(['ENSURE FAILED: ' testCondition]);
            end
        end
        
        function Assert(testCondition)
            return
            if (~evalin('caller', testCondition))
                error(['ASSERT FAILED: ' testCondition]);
            end
        end
        
        function Insist(testCondition, message)
            return
            if (~evalin('caller', testCondition))
                error(['INSIST FAILED: ',testCondition, ' ---> ', message]);
            end
        end
        
    end
    
    methods(Static, Access = private)
        
        function s = join(d, varargin)
        % This JOIN function is an inlined version of Gerald Dalley's one posted at the
        % Matlab Central website.  It is placed here as a convenience to users that
        % have not downloaded it.

            if (isempty(varargin)), 
                s = '';
            else
                if (iscell(varargin{1}))
                    s = DBC.join(d, varargin{1}{:});
                else
                    s = varargin{1};
                end

                for ss = 2:length(varargin)
                    s = [s d DBC.join(d, varargin{ss})];
                end
            end
            
        end
        
        
    end
    
    
end
