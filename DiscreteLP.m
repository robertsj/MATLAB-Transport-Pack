% ======================================================================
%> @brief Compute the DLP's, P_m(i, M) for m = 0 .. M, i = 1 .. M + 1
%
%> m is the "order" and M is the maximum degree.  This follows
%> Neuman, Eqs. 10, 30, and 31.  However, we apply the normalization
%> right away to ease expansions.
% ======================================================================
function [P] = DiscreteLP(M)

    P = zeros(M + 1, M + 1);
    k = (0:M)';         
    
    % Zeroth order
    P(1:end, 1) = 1; 
    
    % First
    if M > 1
        P(1:end, 2) = 1 - (2*k)/M;  
    end
    
    % Higher orders.
    for j = 3:M+1 %otherwise N-1
        i = j-2;
        P(:, i+2) = ( (2*i+1)*(M-2*k).*P(:, i+1) - i*(M+i+1)*P(:, i) ) ./ ...
            ((i+1)*(M-i)) ;
    end
    
    % Normalization
    nrm = sqrt(sum(P.*P));
    P = P*diag(1./nrm);

end