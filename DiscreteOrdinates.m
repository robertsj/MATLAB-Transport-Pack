classdef DiscreteOrdinates
    % Discrete Ordinates solver.
    
    properties (Access = private)
        d_input                 % User input.
        d_mesh                  % Cartesian mesh.
        d_mat                   % Materials.
        d_quadrature            % Angular mesh.
        d_source                % External volume source.
        d_boundary_source       % External boundary source.
        d_fission_source        % Fission source.
        d_state                 % State variables.
        d_solver
    end
    
    methods
        
        function obj = DiscreteOrdinates(input, mesh, mat)
            % function obj = DiscreteOrdinates(input)
            % 
            % Constructor.  Takes an Input class as input.  See Input class for
            % details.
            
            obj.d_input = input;
            obj.d_mesh  = mesh;
            
            % Angular quadrature
            if obj.d_input.quadrature_type == QuadratureTypes.LEVELSYMMETRIC
                obj.d_quadrature = LevelSymmetric(input.quadrature_order);
                
            elseif obj.d_input.quadrature_type == QuadratureTypes.UNIFORMEQUAL
                obj.d_quadrature = UniformEqual(input.quadrature_order);
                
            else
                error('Unknown quadrature type specified.')
                
            end
            
            % State variables
            obj.d_state = State(mesh, input);
            
            
            % Boundary

                                       
            % Fission source is required for all but a fixed source problem
            % without multiplication.  It's only initialized if required.
            %obj.d_fission_source = FissionSource(mesh, mat); 
            
            % Solver
            if strcmp(input.problem_type, 'Fixed')
                obj.d_solver = Fixed(input, mesh, mat, quadrature, boundary);
                
            elseif strcmp(input.problem_type, 'FixedMult')
                obj.d_solver = FixedMult(input, mesh, mat, quadrature);
                
            elseif strcmp(input.problem_type, 'Eigen')
                obj.d_solver = Eigen(input, mesh, mat);
                
            else
               error('Invalid problem_type')
               
            end
            

            
        end
        
        function set_source(obj, source)
        % function set_source(obj, source)
        %   Set an external volume source.
            obj.d_source = source;
        end
        
        function output = solve(obj)
            output = solve(obj.d_solver);
        end
        
    end
    
end


% function [R,F,Ab,LK] = discrete_ordinates(keff, in)
% % function [R,F,Ab,LK] = discrete_ordinates(keff, in)
% %
% % This function solves the 2-D multigroup SN equations given input.  It's 
% % output is the forward (or adjoint) angular fluxes at the cell edges and 
% % the cell-centered scalar flux.  It has been verified against PARTISN for
% % some simple problems.
% %
% % This solver does *not* store angular
% % flux information beyond the last cell in order to reduce the memory load
% % quite substantially.
% %
% % input contains
% %     numg   = number of energy groups
% %     numm   = number of materials
% %     xcm    = x coarse divisions
% %     xfm    = x fine mesh interval for coarse divisions
% %     ycm    = y coarse divisions
% %     yfm    = y fine mesh interval for coarse divisions
% %     mt     = material assignment for each coarse division
% %     data   = cross-sections in the form
% %              (mat1/g1) sigTOT   sigA  sigSg1->g1  sigSg2->g1 ...
% %                 ...g2) sigTOT   sigA  sigSg1->g2  ....
% %              (mat2/g1) ...
% %     src    = volumetric (isotropic) source by coarse mesh and energy
% %              (for adj, src is the cross-section for the source region)
% %     ord    = number of ordinates (2,4,8, or 12)
% %     maxit  = maximum iterations
% %     maxerr = maximum relative pointwise error in phi
% % working variables
% %     dx     = fine mesh divisions
% %     mtt    = material assignment for fine mesh cells
% %     n      = number of fine mesh cells
% 
% numg = in.numg;
% ord  = in.ord;
% [mu,eta,w]=S_2D(ord);
% if ord == 44 || ord == 4.1
%     ord = 4;
% end
% if ord == 66
%     ord = 6;
% end
% if ord == 8.1
%     ord = 8;
% end
% if ord == 16.1
%     ord = 16;
% end
% N = ord;
% if in.quad == 1
%     numord = 0.5*(N+2)*N;
% elseif in.quad == 2
%     numord = N^2;
% else
%     error('Unknown quadrature type!')
% end
% negflx=0;
% 
% % isotropic for now
% Q = zeros(  sum(in.xfm), sum(in.yfm), numg );
% 
% bound       =   [ord/2+1 ord; 1 ord/2];
% gbound      =   [ 1  numg];
% git         =   1;
% 
% if in.adj == 1
%         git         = -1;
%         gbound      = [numg 1];
%         bound       = flipud(bound);
% end
% 
% % -------------------------------------------------------------------------
% % -------------------------------------- Discretizations
% 
% k = 0;
% kk = 0;
% 
% dx = zeros(sum(in.xfm),1);
% dy = zeros(sum(in.yfm),1);
% mtt = zeros(sum(in.xfm),sum(in.yfm));
% idx = zeros(sum(in.xfm),sum(in.yfm),2); % for fine mesh i,j give coarse I,J
% 
% for i = 1:length(in.xfm)
%     dx( (k+1):(k+in.xfm(i))   )  = (in.xcm(i+1)-in.xcm(i))/in.xfm(i);
%     for j = 1:length(in.yfm)
%         dy( (kk+1):(kk+in.yfm(j))   )  = ...
%             (in.ycm(j+1)-in.ycm(j))/in.yfm(j);
%         for g=gbound(1):git:gbound(2)
%             Q( (k+1):(k+in.xfm(i)), (kk+1):(kk+in.yfm(j)), g)  = ...
%                 in.src(g,i,j);
%         end
%         mtt( (k+1):(k+in.xfm(i)), (kk+1):(kk+in.yfm(j))   )  = ...
%             in.mt(i,j);  % assign mat to each f mesh
%         idx( (k+1):(k+in.xfm(i)), (kk+1):(kk+in.yfm(j)),1   )  = i;
%         idx( (k+1):(k+in.xfm(i)), (kk+1):(kk+in.yfm(j)),2   )  = j;       
%         kk = sum(in.yfm(1:j));
%     end
%     kk = 0;
%     k = sum(in.xfm(1:i));
% end
% 
% nx = sum(in.xfm);
% ny = sum(in.yfm);
% 
% % ------------------------------------------------------------------------------
% % Precompute some things.
% 
%   
% % initialize psi matrices (could cut down on storage)
% psiH = zeros(numord,numg);      % horizontal edge flux
% psiV = zeros(ny,numord,numg);   % vertical edge flux
% psiC = zeros(numord,numg);      % cell center
% phi  = zeros(nx,ny,numg); 
% s    = Q;  % isotropic only for now
% 
% % Vacuum boundaries; could easily be updated for reflective
% LeftBoundary    = zeros(ny,numord,numg);
% RightBoundary   = zeros(ny,numord,numg);
% BottomBoundary  = zeros(nx,numord,numg);
% TopBoundary     = zeros(nx,numord,numg);
% 
% % In this memory-efficient implementation, we always go in the vertical
% % direction.  Hence, only the vertical fluxes are stored spatially for use
% % in the next column index
% 
% % ----------------- Solution Algorithm ------------------------------------
% t1=tic;
% for g = gbound(1):git:gbound(2)
%     
%     % ----------------- Convergence Parameters ----------------------------
%     phierr  = 1; % reset phierr
%     it      = 0;
%     A       = 2.0;
%     while it < in.maxit && phierr > in.maxerr
%         phi0(:,:,g) = phi(:,:,g); % store old one
%         phi(:,:,g) = 0; % reset
%         
%         disp([' Starting Iteration: ', num2str(it+1)])
%         % ----------------------------------------------------------------- 
%         % mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
%         % -----------------------------------------------------------------
%         psiV = LeftBoundary; % psiV(:,n,g)
%         for n = 3*numord/4+1:numord
%             for i = 2:nx+1                                  % left to right
%                 psiH(n,g) = BottomBoundary(i-1,n,g); % psiH(n,g)
%                 for j = 2:ny+1                              % bottom to top
%                     m           = mtt(i-1,j-1);
%                     sig         = in.data((m-1)*numg+g,1); % total xsec
%                     con2        = 2*abs(mu(n))/dx(i-1);
%                     con3        = 2*abs(eta(n))/dy(j-1);
%                     con1        = 1 / (sig + con2 + con3);              
%                     psiC(n,g)   = con1 * ( con2*psiV(j-1,n,g) + ...
%                                            con3*psiH(n,g) + ...
%                                            s(i-1,j-1,g) );
%                     if negflx == 1 % simple negative flux fixup              
%                         if psiC(n,g)<0, psiC(n,g)=0; end    
%                     end
%                     psiV(j-1,n,g) = A*psiC(n,g) - psiV(j-1,n,g); % just rewrite over psiV for next column of spatial meshes
%                     psiH(n,g)     = A*psiC(n,g) - psiH(n,g);
%                     if negflx == 1
%                     if psiV(j-1,n,g)<0, psiV(j-1,n,g)=0; end 
%                     if psiH(n,g)<0, psiH(n,g)=0; end 
%                     end
%                     % update the scalar flux
%                     phi(i-1,j-1,g) = phi(i-1,j-1,g) + 0.25*psiC(n,g)*w(n); 
%                     
%                     %fprintf ('%f   %f   %f \n',  ...
%                    % psiC(n,g), psiV(j-1,n,g),psiH(n,g) );                    
%                     
%                 end
%             end
%         end
%         % -----------------------------------------------------------------
%         % mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
%         % -----------------------------------------------------------------
%         psiV = RightBoundary; % psiV(:,n,g)        
%         for n = numord/2+1:3*numord/4
%             for i = nx:-1:1      % right to left
%                 psiH(n,g) = BottomBoundary(i,n,g); % psiH(n,g)                
%                 for j = 2:ny+1   % bottom to top
%                     m           = mtt(i,j-1);
%                     sig         = in.data((m-1)*numg+g,1); % total xsec
%                     con2        = 2*abs(mu(n))/dx(i);
%                     con3        = 2*abs(eta(n))/dy(j-1);
%                     con1        = 1 / (sig + con2 + con3);                     
%                     psiC(n,g) = con1 * ( con2*psiV(j-1,n,g) + ...
%                                          con3*psiH(n,g) + ...
%                                          s(i,j-1,g) );
%                     if negflx == 1
%                     if psiC(n,g)<0, psiC(n,g)=0; end 
%                     end            
%                     psiV(j-1,n,g) = A*psiC(n,g) - psiV(j-1,n,g);
%                     psiH(n,g)     = A*psiC(n,g) - psiH(n,g);
%                     if negflx == 1
%                     if psiV(j-1,n,g)<0, psiV(j-1,n,g)=0; end 
%                     if psiH(n,g)<0, psiH(n,g)=0; end 
%                     end
%                     % update the scalar flux
%                     phi(i,j-1,g) = phi(i,j-1,g) + 0.25*psiC(n,g)*w(n); 
%                     %fprintf ('%f   %f   %f \n',  ...
%                     %psiC(n,g), psiV(j-1,n,g),psiH(n,g) );                    
%                                         
%                 end
%             end
%         end
%         % -----------------------------------------------------------------        
%         % mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
%         % -----------------------------------------------------------------  
%         psiV = LeftBoundary; % psiV(:,n,g)
%         for n = numord/4+1:numord/2
%             for i = 2:nx+1       % left to right
%                 psiH(n,g) = TopBoundary(i-1,n,g); % psiH(n,g)                
%                 for j = ny:-1:1  % top to bottom
%                     m           = mtt(i-1,j);
%                     sig         = in.data((m-1)*numg+g,1); % total xsec
%                     con2        = 2*abs(mu(n))/dx(i-1);
%                     con3        = 2*abs(eta(n))/dy(j);
%                     con1        = 1 / (sig + con2 + con3);                     
%                     psiC(n,g) = con1 * ( con2*psiV(j,n,g) + ...
%                                          con3*psiH(n,g) + ...
%                                          s(i-1,j,g) );
%                     if negflx == 1             
%                     if psiC(n,g)<0, psiC(n,g)=0; end 
%                     end                  
%                     psiV(j,n,g) = A*psiC(n,g) - psiV(j,n,g);
%                     psiH(n,g)   = A*psiC(n,g) - psiH(n,g);
%                     if negflx == 1
%                     if psiV(j,n,g)<0, psiV(j,n,g)=0; end 
%                     if  psiH(n,g)<0, psiH(n,g)=0; end  
%                     end             
%                     % update the scalar flux
%                     phi(i-1,j,g) = phi(i-1,j,g) + 0.25*psiC(n,g)*w(n);       
%                     %fprintf ('%f   %f   %f \n',  ...
%                     %psiC(n,g), psiV(j,n,g),psiH(n,g) );                        
%                 end
%             end
%         end
%         % -----------------------------------------------------------------        
%         % mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
%         % -----------------------------------------------------------------        
%         psiV = RightBoundary; % psiV(:,n,g)        
%         for n = 1:numord/4
%             for i = nx:-1:1      % right to left
%                 psiH(n,g) = TopBoundary(i,n,g); % psiH(n,g)                
%                 for j = ny:-1:1  % top to bottom
%                     m           = mtt(i,j);
%                     sig         = in.data((m-1)*numg+g,1); % total xsec
%                     con2        = 2*abs(mu(n))/dx(i);
%                     con3        = 2*abs(eta(n))/dy(j);
%                     con1        = 1 / (sig + con2 + con3);                    
%                     psiC(n,g) = con1 * ( con2*psiV(j,n,g) + ...
%                                          con3*psiH(n,g) + ...
%                                          s(i,j,g) );   
%                     if negflx == 1           
%                     if psiC(n,g)<0, psiC(n,g)=0; end        
%                     end
%                     psiV(j,n,g) = A*psiC(n,g) - psiV(j,n,g);
%                     psiH(n,g)   = A*psiC(n,g) - psiH(n,g);
%                     if negflx == 1
%                     if psiV(j,n,g)<0, psiV(j,n,g)=0; end 
%                     if psiH(n,g)<0, psiH(n,g)=0; end 
%                     end
%                     % update the scalar flux
%                     phi(i,j,g) = phi(i,j,g) + 0.25*psiC(n,g)*w(n);    
%                     %fprintf ('%f   %f   %f \n',  ...
%                     %psiC(n,g), psiV(j,n,g),psiH(n,g) );                        
%                 end
%             end
%         end    
% %         % ---- Scalar Flux (L&M eq 4-17, l=0)
% %         phi0(:,:) = phi(:,:,g);
% %         phi(:,:,g) = zeros(nx,ny,1);
% %         for n = 1:numord
% %             phi(:,:,g) = phi(:,:,g) + 0.25*psiC(:,:,n,g)*w(n);
% %         end
%         % ---- Updated Source Term
%         for z   =   g:git:gbound(2) % only down scattering
%             for i = 1:nx
%                 for j = 1:ny
%                     s(i,j,z) = Q(i,j,z); % reset
%                 end
%             end
%             if in.adj == 0 % FORWARD TRANSPORT
%                 for i = 1:nx
%                     for j = 1:ny
%                         m = mtt(i,j);
%                         for gg =  gbound(1):git:gbound(2) % group gg to group z
%                            s(i,j,z) = s(i,j,z) + ...
%                                in.data((m-1)*numg+gg,2+z)*phi(i,j,gg);
%                         end
%                     end
%                 end
%             else % ADJOINT TRANSPORT
% %                 for i = 1:nx
% %                     for j = 1:ny
% %                         m = mtt(i,j);
% %                         for gg = gbound(1):git:gbound(2) % group gg to group z
% %                             s(i,j,:,z) = s(i,j,:,z) + ...
% %                                 in.data((m-1)*numg+z,2+gg)*phi(i,j,gg);
% %                         end
% %                     end
% %                 end
%             end
%         end
%         % ---- Scalar Flux Error Between Iterations
%         phierr = max(max(abs(phi(:,:,g)-phi0)./phi(:,:,g)));
%         it = it+1;
%         % Reporting iteration completion
%         fprintf ('Completed inner iteration %g\n',it);
%         fprintf ('Cumulative time %f\n',toc(t1));
%         
%     end % while
%     disp([' Group ', num2str(g),' Iterations:   ',num2str(it)])
% end
% 
% disp([' Elapsed time: ',num2str(toc)])
% disp([' Iterations:   ',num2str(it)])


