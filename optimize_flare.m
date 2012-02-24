%> @file  optimize_flare.m
%> @brief Demonstration optimization of FLARE parameters
% 
%  For core 1
%   -1.1850923847448525   0.03059204859146214   0.27036301829719633
%  For core 2
%   -1.4911192363312082	-0.3564574661838369	0.07210842603953055
function y = optimize_flare(x)
core = 3;
% Get the sample materials for fresh, once-, and twice-burned bundles.
input = Input();
put(input, 'node_width',  21  );
if core == 1 % naive checkerboard
    mat = flare_mat(1);
    % reference values (from a second order response matrix solve)
    refp = [1.8246e+00 2.3824e+00 2.0579e+00 2.0960e+00 1.5013e+00 1.2279e+00 6.1444e-01 1.9250e-01 ...
                       2.1058e+00 2.2998e+00 1.7948e+00 1.6537e+00 1.0342e+00 6.6646e-01 1.8749e-01 ...
                       2.2998e+00 1.8915e+00 1.8804e+00 1.3132e+00 1.0213e+00 4.6492e-01 1.3950e-01 ...
                       1.7948e+00 1.8804e+00 1.4164e+00 1.2478e+00 6.9612e-01 2.6969e-01 8.3160e-02 ...
                       1.6537e+00 1.3132e+00 1.2478e+00 8.2396e-01 5.7262e-01 1.7873e-01 ...
                       1.0342e+00 1.0213e+00 6.9612e-01 5.7262e-01 2.8999e-01 9.2970e-02 ...
                       6.6646e-01 4.6492e-01 2.6969e-01 1.7873e-01 9.2970e-02 ...
                       1.8749e-01 1.3950e-01 8.3160e-02]';
    refmaxpeak = 2.38239;
    refkeff  = 1.179324;
    % The bundle placement.
    map = [ 3  1  2  1  2  1  2  3  0
            0  2  1  2  1  2  1  3  0
            0  1  2  1  2  1  2  3  0
            0  2  1  2  1  2  3  3  0
            0  1  2  1  2  1  3  0  0
            0  2  1  2  1  2  3  0  0
            0  1  2  3  3  3  0  0  0
            0  3  3  3  0  0  0  0  0
            0  0  0  0  0  0  0  0  0];
elseif core == 2 % optimized (mostly; hand tweaking to get symmetry)
    mat = flare_mat(1);
	refp = [ 0.6819   0.87716  1.25513  1.13645  1.50197  1.09531  0.69096  0.30141  ...   
                      0.90524  0.99066  1.45767  1.27607  1.19882  1.08406  0.36121  ...    
                      0.99066  1.36988  1.20131  1.19951  1.52779  1.28003  0.40359  ...    
                      1.45767  1.20131  1.07564  1.08588  1.20188  1.12725  0.33705  ...    
                      1.27607  1.19951  1.08588  1.07886  1.19084  0.81534  ... 
                      1.19882  1.52779  1.20188  1.19084  0.70824  0.30269  ...
                      1.08406  1.28003  1.12725  0.81534  0.30269  ...
                      0.36121  0.40359  0.33705 ]';
    refmaxpeak = 1.52779;
    refkeff    = 1.140545;
    map = [  3     2     1     3     1     2     3     3     0
             0     2     3     1     2     2     1     3     0
             0     3     1     2     2     1     1     3     0
             0     1     2     2     2     2     1     3     0
             0     2     2     2     2     1     1     0     0
             0     2     1     2     1     2     3     0     0
             0     1     1     1     1     3     0     0     0
             0     3     3     3     0     0     0     0     0
             0     0     0     0     0     0     0     0     0];
else
    % 0.9884979114056625	-0.13803858216284481	-0.004050725958076542
    mat = flare_mat(2);
    refp =[ 1.09075   1.10125   1.24251   1.22033   1.08871   0.98202   1.09489   1.01482 ...
                      1.11730   1.13386   1.22347   1.06763   1.03198   1.07188   0.97076 ...
                      1.13386   1.12215   1.10492   1.12001   0.92364   0.93104   0.82465 ...
                      1.22347   1.10492   1.16078   1.03870   0.95013   0.76518   0.54554 ...
                      1.06763   1.12001   1.03870   1.12259   0.99303   0.87444  ...
                      1.03198   0.92364   0.95013   0.99303   1.19960   0.68432   ...
                      1.07188   0.93104   0.76518   0.87444   0.68432   ...
                      0.97076   0.82465   0.54554  ]';
    map = [   1  8  2  6  1  7  1  4  0 
              0  1  8  2  8  1  1  4  0
              0  8  1  8  2  7  1  4  0
              0  2  8  2  8  1  8  4  0
              0  8  2  8  2  5  4  0  0
              0  1  7  1  5  4  4  0  0
              0  1  1  8  4  4  0  0  0
              0  4  4  4  0  0  0  0  0
              0  0  0  0  0  0  0  0  0];
    refkeff = 1.025111;
    refmaxpeak = 1.24251;
    put(input, 'node_width',  23.1226   );
end



% Generate the input and add FLARE-specific things.  

if length(x) == 3
    x(4) = 0;
end
put(input, 'mixing_factor', x(1));
put(input, 'albedo_single', x(2));
put(input, 'albedo_double', x(3)); 
put(input, 'mixing_factor_4', x(4));

% Create the solver.  This generates the coupling coefficients.
solver = FLARE(input, mat, map);

% Solve, returning the fission density and k-eigenvalue.
[s, k] = solve(solver);
k
mean_s = mean(s);
p  = s / mean_s;
max_p  = max(p);

err_k  = abs(k - refkeff)/refkeff
err_p  = norm((p-refp)./refp, Inf)
err_max_p = abs(max_p - refmaxpeak)/refmaxpeak

y = 10*err_k + 100*err_max_p +  10*err_p;

end
