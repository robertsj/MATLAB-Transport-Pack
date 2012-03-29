%> @file  jacobian_matvec.m
%> @brief Jacobian shell matrix-vector operator.
%
%>
%> The Jacobian is defined as 
%>  @f[
%>  \mathbf{f'(x)} = \left [\begin{array}{ccc}
%>	  (\mathbf{M}\mathbf{R}-\lambda \mathbf{I})                       &
%>	            \mathbf{M}\mathbf{R_k}\mathbf{J_-}                     &
%>	                      \mathbf{J_-}                                  \\
%>	  (\mathbf{F}-k\mathbf{L})                                        &
%>	            (\mathbf{F_k}-k\mathbf{L_k}-\mathbf{L}) \mathbf{J_-}   &
%>	                      0                                             \\
%>	  \mathbf{J^T_-}                                                  &
%>	            0                                                      &
%>	                      0
%>	\end{array} 
%>  \right ]  \, .
%>  @f]
%>
%> Note that all but the first @f$ n+1 @f$ rows of the  @f$ n+1 @f$th 
%> column are derivatives with respect to @f$ k @f$. Everything else is 
%> known already.  Hence, the action is computed with a single finite
%> difference approximation.  This is much cheaper than the typical form
%> of JFNK in which each application might involve a different Krylov 
%> vector, the @f$ n+1 @f$th component of which would change and require
%> updated responses.  Here, we need just one update, and beyond that, the
%> Jacobian is analytical.
%>
%> @param   Jac     The shell jacobian object 
%> @param   xx      Incoming PETSc vector
%> @param   yy      Outgoing PETSc vector
%> @param   this    ERME_Newton object making the call
function err = jacobian_matvec(Jac, vv, yy, this)

err = 0;

% Grab incoming vector as MATLAB array.
v = vv(:);
% Decompose for clarity.
v_J       = v(1:end-2, 1);
v_k       = v(end-1);
v_lambda  = v(end);

% Decompose solution
x         = this.d_x(:);
J         = x(1:end-2, 1);
k         = x(end-1);
lambda    = x(end);

n         = length(v);
del_k     = 1e-9;

% Update the operators (won't if k is the same!) and get them.
update(this.problem(), k);
[R1, F1, A1, L1, M, leak] = get_operators(this.problem());

update(this.problem(), k+del_k);
[R2, F2, A2, L2, M, leak] = get_operators(this.problem());

%   [ (M*R - lam*I)    M*R_k*J                                       -J ] [vJ      ]
%Jv=[ F - k*(leak*L+A) F_k*J - leak*L_k*J + A_k*J - (leak*L*J + A*J)  0 ].[v_k     ]
%   [ -J'              0                                              0 ] [v_lambda]

% Compute the action.
y_J         = M*(R1*v_J) - lambda*v_J + ...
              M*((R2*J - R1*J)/del_k)*v_k - ...
              J*v_lambda;
          
y_k         = F1*v_J - k*(leak*(L1*v_J) + A1*v_J) + ...
              ((F2*J-F1*J)/del_k - k*(leak*(L2*J-L1*J)+(A2*J-A1*J))/del_k - ...
                (leak*L1*J + A1*J))*v_k;
              
y_lambda    = -J' * v_J; 


% Fill the outgoing vector
yy(:) = [y_J; y_k; y_lambda];

end