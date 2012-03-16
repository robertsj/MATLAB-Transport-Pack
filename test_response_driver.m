% test of response driver
%clear all
% ==============================================================================
% INPUT 0.836398625920381
% ==============================================================================
input = Input();
put(input, 'number_groups',         1);
put(input, 'dimension',             2);
put(input, 'inner_tolerance',       1e-12);
put(input, 'inner_max_iters',       100);
put(input, 'outer_tolerance',       1e-10);
put(input, 'outer_max_iters',       10);
put(input, 'quad_order',            2);
mso = 0;
mao = 0;
mpo = 0;
for so = mso
for ao = mao
for po = mpo

put(input, 'rf_max_order_space',       so);
put(input, 'rf_max_order_azimuth',     ao);
put(input, 'rf_max_order_polar',       po);
put(input, 'rf_order_space',       so);
put(input, 'rf_order_azimuth',     ao);
put(input, 'rf_order_polar',       po);
kv = linspace(0.5, 1.4, 19);
kv = [ 0.686189327692726]; % 0.0954288045, 0.096151594473819
kv = [ 0.096158350188846]; % flux=0.095536639062328  , bal=0.095521247194524
put(input, 'rf_k_vector',   kv');
put(input, 'rf_number_nodes', 1);
put(input, 'rf_db_name', 'poop.h5');

% ==============================================================================
% MATERIALS (Test two group data)
% ==============================================================================
mat = test_materials(1);

% ==============================================================================
% PINS
% ==============================================================================

% Shared pin properties
pitch = 1.26; radii = 0.54; number = 8;
% Pin 1 - Fuel 1
matid = [2 1];  pin1 = PinCell(pitch, radii, matid); meshify(pin1, number);
% Pin 2 - Fuel 2
matid = [4 1];  pin2 = PinCell(pitch, radii, matid); meshify(pin2, number);
% Pin 3 - MOD
matid = [  1];  pin3 = PinCell(pitch,    [], matid); meshify(pin3, number);
% Pin 4 - Pure Fuel 2
matid = [  2];  pin4 = PinCell(pitch,    [], matid); meshify(pin4, number);

% ==============================================================================
% ASSEMBLIES
% ==============================================================================

% Assembly 1 (kinf =  1.319782)  1.319782 0.891295 1.297046 | 1.278287
pin_map1 = [1 1 1 
            1 1 1
            1 1 1];     
mesh1 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh1);

% Assembly 2 (kinf =  0.891295)
pin_map1 = [1 1 1 
            1 2 1
            1 1 1];     
mesh2 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh2);

% Assembly 3 (kinf =  1.297046)
pin_map1 = [1 1 1 
            1 3 1
            1 1 1];     
mesh3 = Assembly({pin1, pin2, pin3}, pin_map1); 
meshify(mesh3);

% Assembly 4 (kinf =   1.278287)
pin_map1 = [4 4 4
            4 4 4
            4 4 4];     
mesh4 = Assembly({pin1, pin2, pin3, pin4}, pin_map1); 
meshify(mesh4);

% ==============================================================================
% RESPONSE FUNCTION DRIVER
% ==============================================================================
mesh_array = { mesh4};
driver = ResponseDriver(input, mat, mesh_array);
run(driver);


[R, F, A, L] = get_responses(driver);
RR=R(:,:,1,1);
elements = [1 1 1
            1 1 1
            1 1 1];
number_elements = 1;        

tmp = sparse(RR);
Rlist = {tmp;tmp;tmp;tmp;tmp;tmp;tmp;tmp;tmp};
Rblk = blkdiag(Rlist{:});

put(input, 'bc_left',           'vacuum');
put(input, 'bc_right',          'reflect');
put(input, 'bc_bottom',         'reflect');
put(input, 'bc_top',            'reflect');
connect = Connect(input, elements);
[M, Mleak] = build(connect);
% 
% 
% [v,e]=eigs(M*R, 1); J=v; e
% keff = (F'*J)/(A'*J + Mleak*(L'*J))
% eek(so+1, ao+1) = (kv - keff)/kv * 100
% ee(so+1, ao+1) = e;




end
end
end
% 
% for so = 0:mso
%     for ao = 0:mao
%         for po = 0:mpo
%         [R0 F0 A0 L0, M0] = reduce_R(R, F, A, L, M, 1, mso, mao, mpo, so, ao, po);
%         [J, e] = eigs(M0*R0, 1);
%         keff = (F0'*J)/(A0'*J + Mleak*(L0'*J));
%         errk(so+1, ao+1, po+1) = (kv - keff)/kv * 100;
%         lambda(so+1, ao+1, po+1) = e;
%         end
%     end
% end
% errk;
% lambda;
% 
% z = driver.rates.lr / driver.rates.ic(1);
% R2 = [z(1) z(2) z(3) z(3); z(2) z(1) z(3) z(3); z(3) z(3) z(1) z(2); z(3) z(3) z(2) z(1)]; [J2, e2]=eigs(M*R2,1)
% (F'*J2)/(A'*J2 + Mleak*(L'*J2))


[v, e]=eigs(M*R,1); e
J = v(:, 1);
(F'*J)/(A'*J + Mleak*(L'*J))



return

J = sign(J((abs(J)==max(abs(J)))))*J;


FF = F(:, 2);
AA = A(:, 2); 
LL(:, :) = sparse(L(:, :, 2, 1));
Fblk = [FF' FF' FF' FF' FF' FF' FF' FF' FF'];
Ablk = [AA' AA' AA' AA' AA' AA' AA' AA' AA'];
Lblk = {LL' LL' LL' LL' LL' LL' LL' LL' LL'};
Lblk = blkdiag(Lblk{:});



gain        = Fblk*J;                            % compute gains
absorb      = Ablk*J;                            % absorption
leak        = Mleak*(Lblk*J);                     % leakage
loss        = leak + absorb;                  % total loss
gain, absorb, leak, loss, gain/loss

ref.R = Rblk;
ref.M = M;
ref.Mleak = Mleak;
ref.F = Fblk;
ref.A = Ablk;
ref.L = Lblk;
ref.gain = gain;
ref.absorb = absorb;
ref.leak = leak;
ref.loss = loss;
ref.k_input = 0.740382771979217;
ref.k_output = gain/loss;
ref.J = J;

%save('reference.mat', 'ref');

% 
% e = eigs(M*RR, 4, 'LR')
% 
% gain / loss

% % %  
%rf_db = ResponseDB(input);
% initialize_write(rf_db);
% 
% for i = 1:length(mesh_array)
%    write_response(rf_db, i, ['assembly',num2str(i)], ...
%        R(:,:,:,i), F(:,:,i), A(:,:,i), L(:,:,:,i));
% end
% 

% This should yield the SAME
% read_response(rf_db);
% [R, F, A, L] = get_responses(rf_db,  0.740382772061276);
% RR=R(:,:);     
% tmp = sparse(RR);
% Rlist = {tmp;tmp;tmp;tmp;tmp;tmp;tmp;tmp;tmp};
% Rblk = blkdiag(Rlist{:});
% [v, e]=eigs(M*Rblk); 
% J = v(:, 1);
% J = sign(J((abs(J)==max(abs(J)))))*J;
% FF = F;
% AA = A; 
% LL = L;
% Fblk = [FF' FF' FF' FF' FF' FF' FF' FF' FF'];
% Ablk = [AA' AA' AA' AA' AA' AA' AA' AA' AA'];
% Lblk = {LL' LL' LL' LL' LL' LL' LL' LL' LL'};
% Lblk = blkdiag(Lblk{:});
% gain        = Fblk*J;                            % compute gains
% absorb      = Ablk*J;                            % absorption
% leak        = Mleak*(Lblk*J);                     % leakage
% loss        = leak + absorb;                  % total loss
% gain, absorb, leak, loss, gain/loss


