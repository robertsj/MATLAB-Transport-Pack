classdef LevelSymmetric < Quadrature
    % Level Symmetric quadrature.
    %
    % Orders available are 2:2:16.  For higher orders, QuadrupleRange or
    % UniformEqualWeight quadratures are available.
    
    properties
        
    end
    
    methods
        
        % ----------------------------------------------------------------------
        % Constructor
        % ----------------------------------------------------------------------
        
        function obj = LevelSymmetric(order)
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
            % Get data from an old table function.
            [obj.d_mu, obj.d_eta, obj.d_weights] = level_sym_table(order);
            obj.d_number_angles = length(obj.d_weights);
            if (obj.d_number_angles ~= 0.5*(order+2)*order)
                error('Whoops! LevelSymmetric table has an error.')
            end
        end
  
    end
    
    methods (Access = private)
       

    end
    
end

function [mu,eta,w] = level_sym_table(N)
% function [mu,eta,w] = level_sym_table(N)
%   Returns level symmetric quadrature set. 
%  
%   Reference:
%     Lewis and Miller, Computational Methods of Neutron Transport


% note:  N = "44" corresponds to the PARTISN N=4 case, which differs from
% the typical level symmetric LQ values in L&M table 4-1.

% in all cases the order is
%    -mu, -eta
%     mu, -eta
%    -mu,  eta
%     mu,  eta
% moreover, the NN N(N+2)/2 angle pairs (e.g N=4 --> 12 pairs) are 
% arranged such the pair (1+i) is a sign flip of pair (N-i).

if N == 2 
    mus = [ -0.5773503  -0.5773503
             0.5773503  -0.5773503
            -0.5773503   0.5773503
             0.5773503   0.5773503 ];
    w(1:4,1) = 1; 
    
elseif N == 4
    mus = [-0.8688903 -0.3500212   % 1
           -0.3500212 -0.8688903   % 2
           -0.3500212 -0.3500212   % 3
            0.8688903 -0.3500212   % 4
            0.3500212 -0.8688903   % 5
            0.3500212 -0.3500212   % 6
           -0.3500212  0.3500212   % 7
           -0.3500212  0.8688903   % 8
           -0.8688903  0.3500212   % 9    
            0.3500212  0.3500212   %10
            0.3500212  0.8688903   %11
            0.8688903  0.3500212]; %12 
    w(1:12,1) = 1/3;
    
elseif N == 6
    mus = [-0.9261808 -0.2666355 % 1
           -0.6815076 -0.6815076 % 2
           -0.6815076 -0.2666355 % 2
           -0.2666355 -0.9261808 % 1
           -0.2666355 -0.6815076 % 2
           -0.2666355 -0.2666355 % 1
            0.9261808 -0.2666355 % 1
            0.6815076 -0.6815076 % 2
            0.6815076 -0.2666355 % 1
            0.2666355 -0.9261808 % 1
            0.2666355 -0.6815076 % 2
            0.2666355 -0.2666355 % 1
           -0.2666355  0.2666355 % 
           -0.2666355  0.6815076 % 
           -0.2666355  0.9261808 %
           -0.6815076  0.2666355
           -0.6815076  0.6815076  
           -0.9261808  0.2666355
            0.2666355  0.2666355
            0.2666355  0.6815076
            0.2666355  0.9261808
            0.6815076  0.2666355
            0.6815076  0.6815076  
            0.9261808  0.2666355 ];
    w1 = 0.1761263; 
    w2 = 0.1572071;
    w(1:12,1)  = [w1 w2 w2 w1 w2 w1 w1 w2 w1 w1 w2 w1];
    w(13:24,1) = w(12:-1:1,1); 
 
elseif N == 8
    wt1 = 0.1209877; wt2 = 0.0907407; wt3 = 0.0925926;
    mu1 = 0.2182179; mu2 = 0.5773503; mu3 = 0.7867958; mu4 = 0.9511897;
    mu  = [ mu4 mu1    %wt1     411
            mu3 mu2    %wt2     321
            mu3 mu1    %wt2     312
            mu2 mu3    %wt2     231
            mu2 mu2    %wt3     222
            mu2 mu1    %wt2     213
            mu1 mu4    %wt1     141
            mu1 mu3    %wt2     132
            mu1 mu2    %wt2     123
            mu1 mu1 ]; %wt1     114
    mus = zeros(40,2);
    mus(1:10,:)     = -mu;
    mus(11:20,1)    = mu(:,1);  
    mus(11:20,2)    = -mu(:,2);
    mus(21:30,1)    = -mu(:,1); 
    mus(21:30,2)    = mu(:,2);
    mus(31:40,1)    = mu(:,1);  
    mus(31:40,2)    = mu(:,2);
    w(1:10,1)  = [wt1 wt2 wt2 wt2 wt3 wt2 wt1 wt2 wt2 wt1];
    w(11:20,1) = w(1:10,1); 
    w(21:30,1) = w(1:10,1); 
    w(31:40,1) = w(1:10);
    
elseif N == 16
    wt1 = 0.0489872; wt2 = 0.0413296; wt3 = 0.0212326; wt4 = 0.0256207;
    wt5 = 0.0360486; wt6 = 0.0144589; wt7 = 0.0344958; wt8 = 0.0085179;
    mu1 = 0.1389568; mu2 = 0.3922893; mu3 = 0.5370966; mu4 = 0.6504264;
    mu5 = 0.7467506; mu6 = 0.8319966; mu7 = 0.9092855; mu8 = 0.9805009;
    mu = [ mu8 mu1 % 811    wt1
           mu7 mu2 % 721    wt2
           mu7 mu1 % 712    wt2
           mu6 mu3 % 631    wt3
           mu6 mu2 % 622    wt5
           mu6 mu1 % 613    wt3
           mu5 mu4 % 541    wt4
           mu5 mu3 % 532    wt6
           mu5 mu2 % 523    wt6
           mu5 mu1 % 514    wt4
           mu4 mu5 % 451    wt4
           mu4 mu4 % 442    wt7
           mu4 mu3 % 433    wt8
           mu4 mu2 % 424    wt7
           mu4 mu1 % 415    wt4
           mu3 mu6 % 361    wt3
           mu3 mu5 % 352    wt6
           mu3 mu4 % 343    wt8
           mu3 mu3 % 334    wt8
           mu3 mu2 % 325    wt6
           mu3 mu1 % 316    wt3
           mu2 mu7 % 271    wt2
           mu2 mu6 % 262    wt5
           mu2 mu5 % 253    wt6
           mu2 mu4 % 244    wt7
           mu2 mu3 % 235    wt6
           mu2 mu2 % 226    wt5
           mu2 mu1 % 217    wt2
           mu1 mu8 % 181    wt1
           mu1 mu7 % 172    wt2
           mu1 mu6 % 163    wt3
           mu1 mu5 % 154    wt4
           mu1 mu4 % 145    wt4
           mu1 mu3 % 136    wt3
           mu1 mu2 % 127    wt2
           mu1 mu1 % 118    wt1
    ];
    mus = zeros(36*4,2);
    mus(1:36,:)           = -mu;
    mus(1*36+1:2*36,1)    = mu(:,1);  mus(1*36+1:2*36,2)    = -mu(:,2);
    mus(2*36+1:3*36,1)    = -mu(:,1); mus(2*36+1:3*36,2)    = mu(:,2);
    mus(3*36+1:4*36,1)    = mu(:,1);  mus(3*36+1:4*36,2)    = mu(:,2);
    w(1:36,1) = ...
         [ wt1 wt2 wt2 wt3 wt5 wt3 wt4 wt6 wt6 wt4 wt4 wt7 wt8 wt7 wt4 ...
           wt3 wt6 wt8 wt8 wt6 wt3 wt2 wt5 wt6 wt7 wt6 wt5 wt2 wt1 wt2 ...
           wt3 wt4 wt4 wt3 wt2 wt1 ];
    w(1*36+1:2*36,1) = w(1:36,1); 
    w(2*36+1:3*36,1) = w(1:36,1); 
    w(3*36+1:4*36,1) = w(1:36,1);

else
    error('Invalid Level Symmetric order.')
    
end
w = w * pi;    
mu  = -1*mus(:,1); 
eta = -1*mus(:,2);

end
