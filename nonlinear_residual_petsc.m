function err = nonlinear_residual(snes, xx, f, this)
    % Decompose unknown vector for clarity.
    %xx(:)
    x = xx(:);
    J      = x(1:end-2, 1);
    k      = x(end-1);
    lambda = x(end); 

    % Update the operators (won't if k is the same!) and get them.
    update(this.problem(), k);
    [R, F, A, L, M, leak] = get_operators(this.problem());
    % Compute the residual 
    f_J = M*(R*J) - lambda*J;
    f_k = F*J - k*(A*J + leak*(L*J));
    f_lambda = 0.5 - 0.5*(J'*J);
    f(:) = [f_J; f_k; f_lambda];
    err = 0;
end

