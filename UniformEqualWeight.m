classdef UniformEqualWeight < Quadrature
    % Uniform Equal Weight quadrature.
    %
    % Arbitrary order, uniform weight quadrature.  Uniform in azimuth, and
    % Gauss-Legendre in polar.
    %
    % This is actually a product quadrature, since the numbers of polar and
    % azimuthal angles are independent.
    %
    % Reference:
    % 	Carew and Zamonsky, NSE 131 (1999)
    
    properties
        
    end
    
    methods
        
        % ----------------------------------------------------------------------
        % Constructor
        % ----------------------------------------------------------------------
        
        function obj = UniformEqualWeight(order)
            % function obj = Quadrature(order)
            %
            % Constructor.
            %
            % Inputs:
            %   order --    Quadrature order. This differs from quadrature to
            %               quadrature, but e.g. for level symmetric, it's the
            %               number of unique directional cosines.
            obj = obj@Quadrature(order); % Call base class.
            obj.d_order = order;         % Set order.
            % Get data from an old function.  Note, the arguments indicate the
            % total number of angles in 3-D and the number of polar angles,
            % respectively.  Here, we use a descent default that gives a number
            % somewhat close to Level Symmetric.
            [obj.d_mu, obj.d_eta, obj.d_weights] = UEN(2*order^2, order);
            obj.d_number_angles = length(obj.d_weights);
        end
  
    end
    
    methods (Access = private)
       

    end
    
end

function [mu, eta, w] = UEN(N, Nksi)
% function  muwt = UEN(N,Nksi)
% This function computes the Legendre Equal Weights (Pn-EW) Quadrature Set,
% as given by Carew and Zamonsky, NSE 131 (1999).  
%   Inputs:
%       N       -- total number of directions IN 3-D
%       Nksi    -- total number of polar angles IN 3-D
%   Outputs:
%       muwt    -- [mu's eta's wt's] arranged in proper order for code
%
% Note: For 2-d, this could be done a bit more simply by avoiding 
% computation of negative polar and negative.

for i = 1:Nksi % Only compute positive polar
    ksi(i) = (2/Nksi)*(i-1/2)-1;
end

% now for the mu's
delksi = ksi(2)-ksi(1);

Nomega = round(N*delksi/2);
omega = zeros(Nomega, 1);
mu = zeros(Nomega, Nksi);
eta = mu;
for j = 1:Nomega
    omega(j) = (j-0.5)*4*pi/N/delksi;
    for i = 1:Nksi
        mu(i,j)  = cos(omega(j)) * sqrt(1-ksi(i)^2);
        eta(i,j) = sin(omega(j)) * sqrt(1-ksi(i)^2); 
    end
end
% comput vectors of mu's, eta's, weights arranged by
%  [ -mu -eta w
%    +mu -eta w
%    -mu +eta w
%    +mu +eta w ]
k=0;
for j = 1:Nomega/4                  % +mu +eta
    for i = Nksi/2+1:Nksi
        k=k+1;
        muwt(k,:)= [  mu(i,j)  eta(i,j) 8/N ];
    end
end
for j = Nomega/4+1:Nomega/2         % -mu +eta
    for i = Nksi/2+1:Nksi
        k=k+1;
        muwt(k,:)= [ -mu(i,j)  eta(i,j) 8/N ];
    end
end 
for j = Nomega/2+1:3*Nomega/4       % -mu -eta
    for i = Nksi/2+1:Nksi
        k=k+1;
        muwt(k,:)= [ -mu(i,j) -eta(i,j) 8/N ];
    end
end
for j = 3*Nomega/4+1:Nomega         % +mu -eta
    for i = Nksi/2+1:Nksi
        k=k+1;
        muwt(k,:)= [  mu(i,j) -eta(i,j) 8/N ];
    end
end



mu  = muwt(:,1);
eta = muwt(:,2);
w   = muwt(:,3) * pi;

end % function UEN