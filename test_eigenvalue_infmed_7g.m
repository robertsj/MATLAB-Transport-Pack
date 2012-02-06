%> @file  test_eigenvalue_infmed_7g.m
%> @brief Tests the eigensolver for a 7 group infinite medium using 1-d.
% ==============================================================================

clear classes

% ==============================================================================
% INPUT
% ==============================================================================
input = Input();
put(input, 'number_groups',         7);
put(input, 'eigen_tolerance',       1e-4);
put(input, 'eigen_max_iters',       100);
put(input, 'outer_tolerance',       0.0001);
put(input, 'outer_max_iters',       20);
put(input, 'inner_tolerance',       0.0001);
put(input, 'inner_solver',          'SI');
put(input, 'bc_left',               'reflect');
put(input, 'bc_right',              'reflect');
put(input, 'print_out',             1);

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = C5G7_materials(0);

% 
% del.psi(r, Omega) + T*psi(r, Omega) = (S/4pi)*phi + (F/4pi/k)*phi
% Take away spatial dependendence, and integrate over 4pi to get
%   T*phi = S*phi + X*(F/k)*phi
% or
%   inv(T-S)*X*F*phi = k*phi
T = zeros(7);
S = zeros(7);
F = zeros(7);
X = zeros(7);
m = 1;
for g = 1:7
   T(g, g) = sigma_t(mat, m, g); 
   for gp = 1:7 
      % g <-- gp
      S(g, gp) = sigma_s(mat, m, g, gp); 
      % here, need nf(1:ng) in each row so F*phi is [nf'*phi; nf'*phi ...]
      F(gp, g) = nu_sigma_f(mat, m, g);
   end
   X(g, g) = chi(mat, m, g);
end
A = (T-S) \ (X*F);
[v, e] = eig( A );

kinf = max(diag(e));
idx  = find(diag(e)==kinf);
eigv = abs(v(:, idx)  )
disp(['kinf = ',num2str(kinf)])
semilogy(eigv, 'k-o', 'LineWidth', 3)
% grid on, xlabel('group'), ylabel('\phi_g')
%if (get(obj.d_input, 'print_out'))


% ==============================================================================
% MESH 
% ==============================================================================

xcm = [0  1 ];
xfm = [ 200 ];
mt  = [  m  ]; % UO2 
mesh = Mesh1D(xfm, xcm, mt);


% ==============================================================================
% SETUP 
% ==============================================================================
state       = State(input, mesh);
quadrature  = GaussLegendre(64);
boundary    = BoundaryMesh(input, mesh, quadrature);
q_e         = Source(mesh, 2);                  % Not initialized = not used.
q_f         = FissionSource(state, mesh, mat);  % Inititalized = used.
initialize(q_f);

% ==============================================================================
% SOLVE 
% ==============================================================================
solver = Eigensolver(input,         ...
                     state,         ...
                     boundary,      ...
                     mesh,          ...
                     mat,           ...
                     quadrature,	...
                     q_e,           ...
                     q_f);
 
tic
out = solve(solver); 
toc
% 
% % ==============================================================================
% % POSTPROCESS 
% % ==============================================================================
% figure(1)
% clf
% for g = 1:7
%     f = flux(state, g);
%     ff(g,1) = f(1);
%     semilogy(f, 'Color', [rand rand rand], 'LineWidth', 2)
%     grid on
%     hold on
% end
% 
% % This better be constant:
% ff ./ eigv
