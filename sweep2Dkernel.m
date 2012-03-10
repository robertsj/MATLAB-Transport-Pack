%#codegen
function [phi, psi_v, psi_h] = ...
    sweep2Dkernel(phi, psi_v, psi_h, nx, yb, xb, sig, con_x, con_y, s, wt, beta)

% Sweep over all cells.
for j = yb(1):yb(3):yb(2)
    
    psi_v_temp = psi_v(j, :); % psi_v(cells, angles)
    
    for i = xb(1):xb(3):xb(2)
        
        % Set incident angular flux.
        %psi_in = [psi_h(i, :); psi_v_temp];
        
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
        
        % *** Here is where coarse mesh boundary
        %     information could be dealt with.
        
    end
    
    % Save outgoing vertical
    psi_v(j, :) = psi_v_temp;
    
end

end