%> @file  run_tests.m
%> @brief Tests a cell array of evaluation statements.
%> @param tests     Cell array of evaluation statements (strings)
%> @return          Two element array of numbers of tests and failures 
function [nt nf] = run_tests(tests)

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
        disp(['test number: ', num2str(i), ' {', tests{i}, '}',' FAILED.'])
        nf = nf + 1;
    end
end
disp(['*** NUMBER PASSED: ', num2str(nt-nf)])
disp(['*** NUMBER FAILED: ', num2str(nf   )])

end