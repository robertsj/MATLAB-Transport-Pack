%> @file  test_pins.m
%> @brief Construct and return simple pin meshes.
function [assem1, assem2, assem3, assem4] = two_group_benchmark_assemblies()

% Shared pin properties

Mesh(); % @todo Without this, it complains about Mesh.  why?

pin_pitch = 1.26; radii = 0.54; number = 8;
% Pin 1 - Fuel 1
matid = [2 1];  
pin1 = PinCell(pin_pitch, radii, matid); 
meshify(pin1, number);
% Pin 2 - Fuel 2
matid = [4 1];  
pin2 = PinCell(pin_pitch, radii, matid); 
meshify(pin2, number);
% Pin 3 - MOD
matid = [  1];  
pin3 = PinCell(pin_pitch,    [], matid); 
meshify(pin3, number);
% Pin 4 - Pure Fuel 2
matid = [  2];  
pin4 = PinCell(pin_pitch,    [], matid); 
meshify(pin4, number);
pins = {pin1, pin2, pin3, pin4};


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


end

