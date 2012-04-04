%> @file  mr_matvec.m
%> @brief Matrix shell matrix-vector operator for MR.
%
%>
%> @param   MR      The shell for MR
%> @param   xx      Incoming PETSc vector
%> @param   yy      Outgoing PETSc vector
%> @param   this    ERME_Picard object making the call
function err = mr_matvec(Jac, vv, yy, this)

err = 0;

% Grab incoming vector as MATLAB array.
v = vv(:);


this
% Get operators.  This *assumes* they are up-to-date.
[R, F, A, L, M, leak] = get_operators(this.problem());

y = M*(R*v);

% Fill the outgoing vector.
yy(:) = y;


end