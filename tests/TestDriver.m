%> @file  TestDriver.m
%> @brief Routines for doing very simple testing.
classdef TestDriver < handle
    
    properties
        %> Total number of tests performed
        d_number_tests = 0
        %> Total number of failures
        d_number_fails = 0
    end
    
    properties (Constant)
        %> Close enough?
        tolerance = 1e-12; 
    end
    
    methods (Access = public)
        
        %> @brief Tests a cell array of evaluation statements.
        %> @param tests     Cell array of evaluation statements (strings)
        %> @return          Two element array of numbers of tests and failures
        function run_tests(this, tests)
            
            nt = length(tests); % number tests
            nf = 0;             % failure counter
            stack = dbstack(1);
            name = stack.name;  % calling function
            disp(['RUNNING TESTS IN: ', name])
            for i = 1:nt
                try
                    if (~evalin('caller', tests{i}))
                        error(['TEST FAILED: ' tests{i}]);
                    end
                catch
                    disp(['test number: ', num2str(i), ...
                        ' {', tests{i}, '}',' FAILED.'])
                    nf = nf + 1;
                end
            end
            disp(['*** NUMBER PASSED: ', num2str(nt-nf)])
            disp(['*** NUMBER FAILED: ', num2str(nf   )])
            
            % Increment total counters
            this.d_number_tests = this.d_number_tests + nt;
            this.d_number_fails = this.d_number_fails + nf;
        end
        
        function results(this)
            disp('*****************************************')
            disp(['*** TOTAL NUMBER PASSED: ', ...
                  num2str(this.d_number_tests-this.d_number_fails)])
            disp(['*** TOTAL NUMBER FAILED: ', ...
                  num2str(this.d_number_fails )])
            disp('*****************************************')            
        end

    end
    
    methods (Static)
        
        %> @brief  Are two floats close enough?
        function bool = almost_equal(x, y)
            bool = abs(x-y) <= TestDriver.tolerance;
        end
        
    end
end