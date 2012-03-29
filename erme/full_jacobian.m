%> @file  full_jacobian.m
%> @brief Full, possibly approximate, Jacobian matrix for preconditioning.
function Jac = full_jacobian(this)

    % Decompose unknown vector for clarity.
    x       = this.d_x(:);
    J       = x(1:end-2, 1);
    k       = x(end-1);
    lambda  = x(end); 
    n       = length(x);
    del_k   = 1e-9;

    % Update the operators.
    update(this.problem(), k);
    [R1, F1, A1, L1, M, leak] = get_operators(this.problem());

    % Build the approximate Jacobian.
    if this.d_approximate_jacobian == 0
        
        % Include no derivatives with respect to keff.       
        Jac = [ M*R1 - lambda*speye(n-2)	zeros(n-2, 1)        	-J
                F1 - k*(leak*L1+A1)         -(leak*L1*J + A1*J) 	0
                -J'                         0                    	0  ];
    
    else
      
        % Include derivatives with respect to keff.
        update(this.problem(), k+del_k);
        [R2, F2, A2, L2, M, leak] = get_operators(this.problem());

        J_jk = M * ( (R2*J - R1*J)/del_k );
        J_kk = (F2*J-F1*J)/del_k - k*(leak*(L2*J-L1*J) + ...
               (A2*J-A1*J))/del_k - (leak*L1*J + A1*J) ;
        J_kj = F1 - k*(leak*L1+A1);

        Jac = [ M*R1 - lambda*speye(n-2)   J_jk       -J
                J_kj                       J_kk        0
                -J'                        0           0  ];

    end
    
    % There might be other approximations worth considering.
    
end


