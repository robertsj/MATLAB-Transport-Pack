%> @file  test_petsc.m
%> @brief Demo PETSc and SLEPc

% Initialize PETSc and set some typical command line options
%PetscInitialize({'-ksp_monitor','-malloc','-malloc_debug','-malloc_dump'});

% Demo 1
% ------
% 
% Define PETSc Mat and Vec from MATLAB arrays.  First, define a simple 
% finite difference system using some neat built-in MATLAB functions.  (Note,
% the problem is based on the MATLAB 'delsq' demo.)
n = 30;
G = numgrid('L', n);
A_m = delsq(G);
N = sum(G(:) > 0); % Number of interior points
b_m = ones(N, 1);  % Right hand side (Dirichlet condition)
x_m = A_m\b_m;     % Solution via 'backslash'
U = G;
U(G>0) = full(x_m(G(G>0)));
clabel(contour(U));
prism
axis square ij
% 
% % Now do it in PETSc
% 
% %   Create a vector and put values in it
% b = PetscVec();
% b.SetType('seq');
% b.SetSizes(10,10);
% b.SetValues(1:10);
% b.SetValues([1,2],[11.5,12.5],Petsc.ADD_VALUES);
% b.AssemblyBegin();
% b.AssemblyEnd();
% b.View;
% x = b.Duplicate();
% %
% %  Create a matrix and put some values in it
% mat = PetscMat();
% mat.SetType('seqaij');
% mat.SetSizes(10,10,10,10);
% for i=1:10
%   mat.SetValues(i,i,10.0);
% end
% mat.AssemblyBegin(PetscMat.FINAL_ASSEMBLY);
% mat.AssemblyEnd(PetscMat.FINAL_ASSEMBLY);
% mat.View;
% %
% %   Create the linear solver, tell it the matrix to use and solve the system
% ksp = PetscKSP();
% ksp.SetType('gmres');
% ksp.SetOperators(mat,mat,PetscMat.SAME_NONZERO_PATTERN);
% ksp.SetFromOptions();
% ksp.Solve(b,x);
% x.View;
% ksp.View;
% 
% %
% %   Free PETSc objects and shutdown PETSc
% %
% x.Destroy();
% b.Destroy();
% mat.Destroy();
% ksp.Destroy();
% PetscFinalize();
