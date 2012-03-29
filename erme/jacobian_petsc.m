function [flg, err] = jacobian_petsc(snes, xx, JJ, PP, this)
    x = xx(:);
    %disp('jacobian!')
    % Decompose unknown vector for clarity.
    J      = x(1:end-2, 1);
    k      = x(end-1);
    lambda = x(end); 
    n = length(x);
    del_k = 1e-9;

    % Update the operators (won't if k is the same!) and get them.
    update(this.problem(), k);
    [R1, F1, A1, L1, M, leak] = get_operators(this.problem());

    update(this.problem(), k+del_k);
    [R2, F2, A2, L2, M, leak] = get_operators(this.problem());

    %     [  (M*R - lam*I)*J          ] 
    % f = [  F*J - k*(leak*L*J + A*J) ]
    %     [  0.5 - 0.5*J'*J           ]


	  %     [ (M*R - lam*I)       M*R_k*J                                           -J
    % J = [ F - k*(leak*L+A)    F_k*J - leak*L_k*J + A_k*J - (leak*L*J + A*J)      0
    %     [ -J'                 0                                                  0

    
    J_jk = M * ( (R2*J - R1*J)/del_k );
    J_kk = (F2*J-F1*J)/del_k - k*(leak*(L2*J-L1*J) + (A2*J-A1*J))/del_k - (leak*L1*J + A1*J) ;
    J_kj = F1 - k*(leak*L1+A1);
 
    Jac = [ M*R1 - lambda*speye(n-2)   J_jk       -J
            J_kj                       J_kk        0
            -J'                        0           0  ];

    save('jac.mat', 'Jac')
    for j = 1:n
      for i = 1:n
        if (Jac(i,j)~=0)
          JJ.SetValues(i,j,Jac(i,j));
        end
      end
    end
    JJ.SetValues(n,n,0.0);
    %diag(Jac), size(Jac)
    %if this.d_built == 0
    %JJ = PetscMat(Jac);
    %PP = PetscMat(speye(n));
    %this.d_built = 1;
    %else
    %disp('something new!')
    %JJ(:, :) = Jac;
    %PP(:, :) = speye(n);
    %end

    err = JJ.AssemblyBegin(PetscMat.FINAL_ASSEMBLY);
    err = JJ.AssemblyEnd(PetscMat.FINAL_ASSEMBLY);
    %err = PP.AssemblyBegin(PetscMat.FINAL_ASSEMBLY);
    %err = PP.AssemblyEnd(PetscMat.FINAL_ASSEMBLY);

    flg = PetscMat.SAME_NONZERO_PATTERN;
    %disp('jac lala')
end

%function [flg,err] = nljacobian(snes,x,A,B,ctx)
%%
%%  Example of a nonlinear Jacobian needed by SNES
%%  Use identity as approximation for Jacobian
%%
%disp('jac')
%disp(['ctx.a=',num2str(ctx.a)])
%err = 0;
%flg = PetscMat.SAME_NONZERO_PATTERN;
%for i=1:10
%  B.SetValues(i,i,2.0);
%end
%err = B.AssemblyBegin(PetscMat.FINAL_ASSEMBLY);
%err = B.AssemblyEnd(PetscMat.FINAL_ASSEMBLY);
%err = A.AssemblyBegin(PetscMat.FINAL_ASSEMBLY);
%err = A.AssemblyEnd(PetscMat.FINAL_ASSEMBLY);

