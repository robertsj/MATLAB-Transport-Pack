%> @file  tests/regression/test_assemblies.m
%> @brief Construct and return simple pin meshes.
function assemblies = test_assemblies(tester, pins)

% Shared pin properties

Mesh(); % @todo Without this, it complains about Mesh.  why?

% Assemblies
pin_map1 = [1 1 1 ; 1 1 1 ; 1 1 1]; % 3x3 UO2 pins
assem1   = Assembly(pins, pin_map1); 
meshify(assem1);
pin_map2 = [1 1 1 ; 1 2 1 ; 1 1 1]; % 3x3 UO2 with Gd type I central
assem2   = Assembly(pins, pin_map2); 
meshify(assem2);
pin_map3 = [1 1 1 ; 1 3 1 ; 1 1 1]; % 3x3 UO2 with Gd type II central
assem3   = Assembly(pins, pin_map3); 
meshify(assem3);
pin_map4 = [4 4 4 ; 4 4 4 ; 4 4 4]; % Homogeneous Gd type I fuel
assem4   = Assembly(pins, pin_map4); 
meshify(assem4);
assemblies = {assem1, assem2, assem3, assem4};

% Define evaluation statements and test tem.
tests = {'tester.almost_equal(pitch(assem4), 3.78)'};
tester.run_tests(tests);

end

