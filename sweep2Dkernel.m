%> @file  sweep2Dkernel.m
%> @brief Performs a 2-D sweep over all angles within an octant.
%
%> This function is intended for use with MEX.  Initial studies indicate
%> it provides a speedup of 3-5 over the pure MATLAB version.
%#codegen
function [phi, psi_v, psi_h] = ...
    sweep2Dkernel(phi, psi_v, psi_h, nx, yb, xb, sig, con_x, con_y, s, wt, beta)

% Sweep over all cells.
for j = yb(1):yb(3):yb(2)
    
    psi_v_temp = psi_v(j, :); % psi_v(cells, angles)
    
    for i = xb(1):xb(3):xb(2)
        
        % Cardinal index.
        k = i + (j - 1) * nx;
        
        coef = 1.0 ./ (sig(i, j) + con_x(i, :) + con_y(j, :));
        psi_center = coef .* ...
            (s(k) + con_x(i, :) .* psi_v_temp + ...
            con_y(j, :) .* psi_h(i, :) );
        
        % Outgoing fluxes
        psi_h(i, :) = beta(1)*psi_center + beta(2)*psi_h(i, :);
        psi_v_temp  = beta(1)*psi_center + beta(2)*psi_v_temp;
        
        % Inner product of weights with psi.
        phi(k) = phi(k) + psi_center * wt;
        
    end
    
    % Save outgoing vertical
    psi_v(j, :) = psi_v_temp;
    
end

end