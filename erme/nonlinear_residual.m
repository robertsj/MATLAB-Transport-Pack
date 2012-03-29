%> @file  nonlinear_residual.m
%> @brief Nonlinear residual for ERME.
%
%> The nonlinear residual is defined as
%> @f[
%>       \mathbf{f(x)} = \left [\begin{array}{c}
%>	        (\mathbf{M}\mathbf{R}(k)-\lambda \mathbf{I}) \mathbf{J_-} \\
%>	        \mathbf{F}(k)\mathbf{J_-} - (k\mathbf{L}(k)\mathbf{J_-} ) \\
%>	        \frac{1}{2} \mathbf{J^T_-} \mathbf{J_-} - \frac{1}{2}  
%>	      \end{array} 
%>       \right ]  = \mathbf{0} \, .
%> @f]
%> 
%> @param   snes    The PETSc SNES object.
%> @param   xx      Incoming PETSc vector
%> @param   f       Outgoing PETSc vector
%> @param   this    ERME_Newton object making the call
function err = nonlinear_residual(snes, xx, f, this)

    % Decompose unknown vector.
    x       = xx(:);
    J       = x(1:end-2, 1);
    k       = x(end-1);
    lambda  = x(end); 

    % Update and retrieve the operators.
    update(this.problem(), k);
    [R, F, A, L, M, leak] = get_operators(this.problem());
    
    % Compute the residual 
    f_J         = M*(R*J) - lambda*J;
    f_k         = F*J - k*(A*J + leak*(L*J));
    f_lambda    = 0.5 - 0.5*(J'*J);
    f(:)      	= [f_J; f_k; f_lambda];
    err         = 0;
end

