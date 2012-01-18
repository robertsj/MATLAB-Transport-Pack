function mat = test_materials(n)

if nargin == 0
    n = 1;
end

% ------------------------------------------------------------------------------
% Simple 1 group, 2 material
% ------------------------------------------------------------------------------
if n == 1
    
% Create the Materials object.
mat = Materials(1, ...  % groups
                2);     % materials         
% Total
mat = set_sigma_t(mat, 1, 1, 1.0); % (object, material, group, value)
% Fission 
mat = set_nu_sigma_f(mat, 1, 1, 0.0); 
mat = set_chi(mat, 1, 1, 0.0);     
% Scattering
mat = set_sigma_s(mat, 1, 1, 1, 0.8); % (object, material, g, g', value);

% Total
mat = set_sigma_t(mat, 2, 1, 1.0); % (object, material, group, value)
% Fission 
mat = set_nu_sigma_f(mat, 2, 1, 0.0); 
mat = set_chi(mat, 2, 1, 0.0);     
% Scattering
mat = set_sigma_s(mat, 2, 1, 1, 0.0); % (object, material, g, g', value);

% ------------------------------------------------------------------------------
elseif n == 2
    
% Create the Materials object.
mat = Materials(2, ...  % groups
                4);     % materials

% ---------------------------
% Material 1: Water           
% ---------------------------

% Total
mat = set_sigma_t(mat, 1, 1, 0.1890); % (object, material, group, value)
mat = set_sigma_t(mat, 1, 2, 1.4633);       

% Fission 
mat = set_nu_sigma_f(mat, 1, 1, 0.0); % Note, default is zero
mat = set_nu_sigma_f(mat, 1, 2, 0.0);   
mat = set_chi(mat, 1, 1, 0.0); 
mat = set_chi(mat, 1, 2, 0.0);        

% Scattering
mat = set_sigma_s(mat, 1, 1, 1, 0.1507); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 1, 1, 2, 0.0000); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 1, 2, 1, 0.0380); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 1, 2, 2, 1.4536); % (object, material, g, g', value);g<-g'.

% ---------------------------
% Material 2: Fuel I           
% ---------------------------

% Total
mat = set_sigma_t(mat, 2, 1, 0.2263);   % (object, material, group, value)
mat = set_sigma_t(mat, 2, 2, 1.0119);       

% Fission 
mat = set_nu_sigma_f(mat, 2, 1, 0.0067);
mat = set_nu_sigma_f(mat, 2, 2, 0.1241);   
mat = set_chi(mat, 2, 1, 1.0); 
mat = set_chi(mat, 2, 2, 0.0);        

% Scattering
mat = set_sigma_s(mat, 2, 1, 1, 0.2006); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 2, 1, 2, 0.0);    % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 2, 2, 1, 0.0161); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 2, 2, 2, 0.9355); % (object, material, g, g', value);g<-g'.

% ---------------------------
% Material 3: Fuel II          
% ---------------------------
%                 0.2252 0.0000 0.0078 1.0000 0.1995 0.0156 % fuel II
%                 0.9915 0.0000 0.1542 0.0000 0.0000 0.9014

% Total
mat = set_sigma_t(mat, 3, 1, 0.2252);   % (object, material, group, value)
mat = set_sigma_t(mat, 3, 2, 0.9915);       

% Fission 
mat = set_nu_sigma_f(mat, 3, 1, 0.0078);
mat = set_nu_sigma_f(mat, 3, 2, 0.1542);   
mat = set_chi(mat, 3, 1, 1.0); 
mat = set_chi(mat, 3, 2, 0.0);        

% Scattering
mat = set_sigma_s(mat, 3, 1, 1, 0.1995); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 3, 1, 2, 0.0000); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 3, 2, 1, 0.0156); % (object, material, g, g', value);g<-g'.
mat = set_sigma_s(mat, 3, 2, 2, 0.9014); % (object, material, g, g', value);g<-g'.
       
% ---------------------------
% Material 4: Fuel III           
% ---------------------------
%                 0.2173 0.0000 0.0056 1.0000 0.1902 0.0136 % fuel III
%                 1.0606 0.0000 0.0187 0.0000 0.0000 0.5733

% Total
mat = set_sigma_t(mat, 4, 1, 0.2173);   % (object, material, group, value)
mat = set_sigma_t(mat, 4, 2, 1.0606);       

% Fission 
mat = set_nu_sigma_f(mat, 4, 1, 0.0056);
mat = set_nu_sigma_f(mat, 4, 2, 0.0187);   
mat = set_chi(mat, 4, 1, 1.0); 
mat = set_chi(mat, 4, 2, 0.0);        

% Scattering
mat = set_sigma_s(mat, 4, 1, 1, 0.1902); % (object, mat, g, g', val);g<-g'.
mat = set_sigma_s(mat, 4, 1, 2, 0.0000); % (object, mat, g, g', val);g<-g'.
mat = set_sigma_s(mat, 4, 2, 1, 0.0136); % (object, mat, g, g', val);g<-g'.
mat = set_sigma_s(mat, 4, 2, 2, 0.5733); % (object, mat, g, g', val);g<-g'.

% ---------------------------
% FINALIZE     
% ---------------------------

% This sets the scattering bounds, which can eliminate a few operations.
mat = finalize(mat);

% ---------------------------
% VERIFY          
% ---------------------------
%DBC.Require('mat.sigma_t(1, 1) == 0.1891');
%DBC.Insist('mat.sigma_t(1, 1) == 0.1891', 'I insist!');

if (mat.sigma_t(1, 1) ~= 0.1890)
    error('Wrong!')
end
if (mat.nu_sigma_f(2, 2) ~= 0.1241)
    error('Wrong!')
end
if (mat.chi(3, 2) ~= 0.0)
    error('Wrong!')
end
if (mat.sigma_s(4, 1, 2) ~= 0)
    error('Wrong!')
end
if (mat.sigma_s(4, 2, 1) ~= 0.0136)
    error('Wrong!')
end

end


end