function mat = C5G7_materials()

% Create the Materials object.
mat = Materials(7, ...  % groups
                7);     % materials

% -----------------------------------------------------------------------------
% (a) UO2 fuel-clad
% -----------------------------------------------------------------------------
m = 1;
% Transport cross section
set_sigma_t(mat, m, 1, 1.77949E-01);
set_sigma_t(mat, m, 2, 3.29805E-01);
set_sigma_t(mat, m, 3, 4.80388E-01);
set_sigma_t(mat, m, 4, 5.54367E-01);
set_sigma_t(mat, m, 5, 3.11801E-01);
set_sigma_t(mat, m, 6, 3.95168E-01);
set_sigma_t(mat, m, 7, 5.64406E-01);
% Absorption cross section
% set_sigma_a(mat, m, 1, 8.02480E-03);
% set_sigma_a(mat, m, 2, 3.71740E-03);
% set_sigma_a(mat, m, 3, 2.67690E-02);
% set_sigma_a(mat, m, 4, 9.62360E-02);
% set_sigma_a(mat, m, 5, 3.00200E-02);
% set_sigma_a(mat, m, 6, 1.11260E-01);
% set_sigma_a(mat, m, 7, 2.82780E-01);
% Fission times nu
set_nu_sigma_f(mat, m, 1, 7.21206E-03*2.78145E+00);
set_nu_sigma_f(mat, m, 1, 8.19301E-04*2.47443E+00);
set_nu_sigma_f(mat, m, 1, 6.45320E-03*2.43383E+00);
set_nu_sigma_f(mat, m, 1, 1.85648E-02*2.43380E+00);
set_nu_sigma_f(mat, m, 1, 1.78084E-02*2.43380E+00);
set_nu_sigma_f(mat, m, 1, 8.30348E-02*2.43380E+00);
set_nu_sigma_f(mat, m, 1, 2.16004E-01*2.43380E+00);
% Fission spectrum
set_chi(mat, m, 1, 5.87819E-01);
set_chi(mat, m, 2, 4.11760E-01);
set_chi(mat, m, 3, 3.39060E-04);
set_chi(mat, m, 4, 1.17610E-07);
set_chi(mat, m, 5, 0.00000E+00);
set_chi(mat, m, 6, 0.00000E+00);
set_chi(mat, m, 7, 0.00000E+00);
% Scattering
%   1 <- g'
set_sigma_s(mat, m, 1, 1, 1.27537E-01);
%   2 <- g'
set_sigma_s(mat, m, 2, 1, 4.23780E-02);
set_sigma_s(mat, m, 2, 2, 3.24456E-01); 
%   3 <- g'
set_sigma_s(mat, m, 3, 1, 9.43740E-06);
set_sigma_s(mat, m, 3, 2, 1.63140E-03);
set_sigma_s(mat, m, 3, 3, 4.50940E-01); 
%   4 <- g'
set_sigma_s(mat, m, 4, 1, 5.51630E-09);
set_sigma_s(mat, m, 4, 2, 3.14270E-09);
set_sigma_s(mat, m, 4, 3, 2.67920E-03); 
set_sigma_s(mat, m, 4, 4, 4.52565E-01);  
set_sigma_s(mat, m, 4, 5, 1.25250E-04);  
%   5 <- g'
set_sigma_s(mat, m, 5, 4, 5.56640E-03);  
set_sigma_s(mat, m, 5, 5, 2.71401E-01);  
set_sigma_s(mat, m, 5, 6, 1.29680E-03);   
%   6 <- g'
set_sigma_s(mat, m, 6, 5, 1.02550E-02);  
set_sigma_s(mat, m, 6, 6, 2.65802E-01);   
set_sigma_s(mat, m, 6, 7, 8.54580E-03);  
%   7 <- g'
set_sigma_s(mat, m, 7, 5, 1.00210E-08);  
set_sigma_s(mat, m, 7, 6, 1.68090E-02);   
set_sigma_s(mat, m, 7, 7, 2.73080E-01);  

% ---------------------------
% (b) 4.3% MOX fuel-clad
% ---------------------------
m = 2;
% Transport cross section
set_sigma_t(mat, m, 1, 1.78731E-01);
set_sigma_t(mat, m, 2, 3.30849E-01);
set_sigma_t(mat, m, 3, 4.83772E-01);
set_sigma_t(mat, m, 4, 5.66922E-01);
set_sigma_t(mat, m, 5, 4.26227E-01);
set_sigma_t(mat, m, 6, 6.78997E-01);
set_sigma_t(mat, m, 7, 6.82852E-01);
% Absorption cross section
% set_sigma_a(mat, m, 1, 8.43390E-03);
% set_sigma_a(mat, m, 2, 3.75770E-03);
% set_sigma_a(mat, m, 3, 2.79700E-02);
% set_sigma_a(mat, m, 4, 1.04210E-01);
% set_sigma_a(mat, m, 5, 1.39940E-01);
% set_sigma_a(mat, m, 6, 4.09180E-01);
% set_sigma_a(mat, m, 7, 4.09350E-01);
% Fission times nu
set_nu_sigma_f(mat, m, 1, 7.62704E-03*2.85209E+00);
set_nu_sigma_f(mat, m, 1, 8.76898E-04*2.89099E+00);
set_nu_sigma_f(mat, m, 1, 5.69835E-03*2.85486E+00);
set_nu_sigma_f(mat, m, 1, 2.28872E-02*2.86073E+00);
set_nu_sigma_f(mat, m, 1, 1.07635E-02*2.85447E+00);
set_nu_sigma_f(mat, m, 1, 2.32757E-01*2.86415E+00);
set_nu_sigma_f(mat, m, 1, 2.48968E-01*2.86780E+00);
% Fission spectrum
set_chi(mat, m, 1, 5.87819E-01);
set_chi(mat, m, 2, 4.11760E-01);
set_chi(mat, m, 3, 3.39060E-04);
set_chi(mat, m, 4, 1.17610E-07);
set_chi(mat, m, 5, 0.00000E+00);
set_chi(mat, m, 6, 0.00000E+00);
set_chi(mat, m, 7, 0.00000E+00);
% Scattering
%   1 <- g'
set_sigma_s(mat, m, 1, 1, 1.28876E-01);
%   2 <- g'
set_sigma_s(mat, m, 2, 1, 4.14130E-02);
set_sigma_s(mat, m, 2, 2, 3.25452E-01); 
%   3 <- g'
set_sigma_s(mat, m, 3, 1, 8.22900E-06);
set_sigma_s(mat, m, 3, 2, 1.63950E-03);
set_sigma_s(mat, m, 3, 3, 4.53188E-01); 
%   4 <- g'
set_sigma_s(mat, m, 4, 1, 5.04050E-09);
set_sigma_s(mat, m, 4, 2, 1.59820E-09);
set_sigma_s(mat, m, 4, 3, 2.61420E-03); 
set_sigma_s(mat, m, 4, 4, 4.57173E-01);  
set_sigma_s(mat, m, 4, 5, 1.60460E-04);  
%   5 <- g'
set_sigma_s(mat, m, 5, 4, 5.53940E-03);  
set_sigma_s(mat, m, 5, 5, 2.76814E-01);  
set_sigma_s(mat, m, 5, 6, 2.00510E-03);   
%   6 <- g'
set_sigma_s(mat, m, 6, 5, 9.31270E-03);  
set_sigma_s(mat, m, 6, 6, 2.52962E-01);   
set_sigma_s(mat, m, 6, 7, 8.49480E-03);  
%   7 <- g'
set_sigma_s(mat, m, 7, 5, 9.16560E-09);  
set_sigma_s(mat, m, 7, 6, 1.48500E-02);   
set_sigma_s(mat, m, 7, 7, 2.65007E-01);  

% ---------------------------
% (c) 7.0% MOX fuel-clad
% ---------------------------
m = 3;
% Transport cross section
set_sigma_t(mat, m, 1, 1.81323E-01);
set_sigma_t(mat, m, 2, 3.34368E-01);
set_sigma_t(mat, m, 3, 4.93785E-01);
set_sigma_t(mat, m, 4, 5.91216E-01);
set_sigma_t(mat, m, 5, 4.74198E-01);
set_sigma_t(mat, m, 6, 8.33601E-01);
set_sigma_t(mat, m, 7, 8.53603E-01);
% Absorption cross section
% set_sigma_a(mat, m, 1, 9.06570E-03);
% set_sigma_a(mat, m, 2, 4.29670E-03);
% set_sigma_a(mat, m, 3, 3.28810E-02);
% set_sigma_a(mat, m, 4, 1.22030E-01);
% set_sigma_a(mat, m, 5, 1.82980E-01);
% set_sigma_a(mat, m, 6, 5.68460E-01);
% set_sigma_a(mat, m, 7, 5.85210E-01);
% Fission times nu
set_nu_sigma_f(mat, m, 1, 8.25446E-03*2.88498E+00);
set_nu_sigma_f(mat, m, 1, 1.32565E-03*2.91079E+00);
set_nu_sigma_f(mat, m, 1, 8.42156E-03*2.86574E+00);
set_nu_sigma_f(mat, m, 1, 3.28730E-02*2.87063E+00);
set_nu_sigma_f(mat, m, 1, 1.59636E-02*2.86714E+00);
set_nu_sigma_f(mat, m, 1, 3.23794E-01*2.86658E+00);
set_nu_sigma_f(mat, m, 1, 3.62803E-01*2.87539E+00);
% Fission spectrum
set_chi(mat, m, 1, 5.87819E-01);
set_chi(mat, m, 2, 4.11760E-01);
set_chi(mat, m, 3, 3.39060E-04);
set_chi(mat, m, 4, 1.17610E-07);
set_chi(mat, m, 5, 0.00000E+00);
set_chi(mat, m, 6, 0.00000E+00);
set_chi(mat, m, 7, 0.00000E+00);
% Scattering
%   1 <- g'
set_sigma_s(mat, m, 1, 1, 1.30457E-01);
%   2 <- g'
set_sigma_s(mat, m, 2, 1, 4.17920E-02);
set_sigma_s(mat, m, 2, 2, 3.28428E-01); 
%   3 <- g'
set_sigma_s(mat, m, 3, 1, 8.51050E-06);
set_sigma_s(mat, m, 3, 2, 1.64360E-03);
set_sigma_s(mat, m, 3, 3, 4.58371E-01); 
%   4 <- g'
set_sigma_s(mat, m, 4, 1, 5.13290E-09);
set_sigma_s(mat, m, 4, 2, 2.20170E-09);
set_sigma_s(mat, m, 4, 3, 2.53310E-03); 
set_sigma_s(mat, m, 4, 4, 4.63709E-01);  
set_sigma_s(mat, m, 4, 5, 1.76190E-04);  
%   5 <- g'
set_sigma_s(mat, m, 5, 4, 5.47660E-03);  
set_sigma_s(mat, m, 5, 5, 2.82313E-01);  
set_sigma_s(mat, m, 5, 6, 2.27600E-03);   
%   6 <- g'
set_sigma_s(mat, m, 6, 5, 8.72890E-03);  
set_sigma_s(mat, m, 6, 6, 2.49751E-01);   
set_sigma_s(mat, m, 6, 7, 8.86450E-03);  
%   7 <- g'
set_sigma_s(mat, 3, 7, 5, 9.00160E-09);  
set_sigma_s(mat, 3, 7, 6, 1.31140E-02);   
set_sigma_s(mat, 3, 7, 7, 2.59529E-01);  

% ---------------------------
% (d) 8.7% MOX fuel-clad
% ---------------------------
m = 4;
% Transport cross section
set_sigma_t(mat, m, 1, 1.83045E-01);
set_sigma_t(mat, m, 2, 3.36705E-01);
set_sigma_t(mat, m, 3, 5.00507E-01);
set_sigma_t(mat, m, 4, 6.06174E-01);
set_sigma_t(mat, m, 5, 5.02754E-01);
set_sigma_t(mat, m, 6, 9.21028E-01);
set_sigma_t(mat, m, 7, 9.55231E-01);
% Absorption cross section
% set_sigma_a(mat, m, 1, 9.48620E-03);
% set_sigma_a(mat, m, 2, 4.65560E-03);
% set_sigma_a(mat, m, 3, 3.62400E-02);
% set_sigma_a(mat, m, 4, 1.32720E-01);
% set_sigma_a(mat, m, 5, 2.08400E-01);
% set_sigma_a(mat, m, 6, 6.58700E-01);
% set_sigma_a(mat, m, 7, 6.90170E-01);
% Fission times nu
set_nu_sigma_f(mat, m, 1, 8.67209E-03*2.90426E+00);
set_nu_sigma_f(mat, m, 1, 1.62426E-03*2.91795E+00);
set_nu_sigma_f(mat, m, 1, 1.02716E-02*2.86986E+00);
set_nu_sigma_f(mat, m, 1, 3.90447E-02*2.87491E+00);
set_nu_sigma_f(mat, m, 1, 1.92576E-02*2.87175E+00);
set_nu_sigma_f(mat, m, 1, 3.74888E-01*2.86752E+00);
set_nu_sigma_f(mat, m, 1, 4.30599E-01*2.87808E+00);
% Fission spectrum
set_chi(mat, m, 1, 5.87819E-01);
set_chi(mat, m, 2, 4.11760E-01);
set_chi(mat, m, 3, 3.39060E-04);
set_chi(mat, m, 4, 1.17610E-07);
set_chi(mat, m, 5, 0.00000E+00);
set_chi(mat, m, 6, 0.00000E+00);
set_chi(mat, m, 7, 0.00000E+00);
% Scattering
%   1 <- g'
set_sigma_s(mat, m, 1, 1, 1.31504E-01);
%   2 <- g'
set_sigma_s(mat, m, 2, 1, 4.20460E-02);
set_sigma_s(mat, m, 2, 2, 3.30403E-01); 
%   3 <- g'
set_sigma_s(mat, m, 3, 1, 8.69720E-06);
set_sigma_s(mat, m, 3, 2, 1.64630E-03);
set_sigma_s(mat, m, 3, 3, 4.61792E-01); 
%   4 <- g'
set_sigma_s(mat, m, 4, 1, 5.19380E-09);
set_sigma_s(mat, m, 4, 2, 2.60060E-09);
set_sigma_s(mat, m, 4, 3, 2.47490E-03); 
set_sigma_s(mat, m, 4, 4, 4.68021E-01);  
set_sigma_s(mat, m, 4, 5, 1.85970E-04);  
%   5 <- g'
set_sigma_s(mat, m, 5, 4, 5.43300E-03);  
set_sigma_s(mat, m, 5, 5, 2.85771E-01);  
set_sigma_s(mat, m, 5, 6, 2.39160E-03);   
%   6 <- g'
set_sigma_s(mat, m, 6, 5, 8.39730E-03);  
set_sigma_s(mat, m, 6, 6, 2.47614E-01);   
set_sigma_s(mat, m, 6, 7, 8.96810E-03);  
%   7 <- g'
set_sigma_s(mat, m, 7, 5, 8.92800E-09);  
set_sigma_s(mat, m, 7, 6, 1.23220E-02);   
set_sigma_s(mat, m, 7, 7, 2.56093E-01);  

% ---------------------------
% (e) fission chamber
% ---------------------------
m = 5;
% Transport cross section
set_sigma_t(mat, m, 1, 1.26032E-01);
set_sigma_t(mat, m, 2, 2.93160E-01);
set_sigma_t(mat, m, 3, 2.84250E-01);
set_sigma_t(mat, m, 4, 2.81020E-01);
set_sigma_t(mat, m, 5, 3.34460E-01);
set_sigma_t(mat, m, 6, 5.65640E-01);
set_sigma_t(mat, m, 7, 1.17214E+00);
% Absorption cross section
% set_sigma_a(mat, m, 1, 5.11320E-04);
% set_sigma_a(mat, m, 2, 7.58130E-05);
% set_sigma_a(mat, m, 3, 3.16430E-04);
% set_sigma_a(mat, m, 4, 1.16750E-03);
% set_sigma_a(mat, m, 5, 3.39770E-03);
% set_sigma_a(mat, m, 6, 9.18860E-03);
% set_sigma_a(mat, m, 7, 2.32440E-02);
% Fission times nu
set_nu_sigma_f(mat, m, 1, 4.79002E-09*2.76283E+00);
set_nu_sigma_f(mat, m, 1, 5.82564E-09*2.46239E+00);
set_nu_sigma_f(mat, m, 1, 4.63719E-07*2.43380E+00);
set_nu_sigma_f(mat, m, 1, 5.24406E-06*2.43380E+00);
set_nu_sigma_f(mat, m, 1, 1.45390E-07*2.43380E+00);
set_nu_sigma_f(mat, m, 1, 7.14972E-07*2.43380E+00);
set_nu_sigma_f(mat, m, 1, 2.08041E-06*2.43380E+00);
% Fission spectrum
set_chi(mat, m, 1, 5.87819E-01);
set_chi(mat, m, 2, 4.11760E-01);
set_chi(mat, m, 3, 3.39060E-04);
set_chi(mat, m, 4, 1.17610E-07);
set_chi(mat, m, 5, 0.00000E+00);
set_chi(mat, m, 6, 0.00000E+00);
set_chi(mat, m, 7, 0.00000E+00);
% Scattering
%   1 <- g'
set_sigma_s(mat, m, 1, 1, 6.61659E-02);
%   2 <- g'
set_sigma_s(mat, m, 2, 1, 5.90700E-02);
set_sigma_s(mat, m, 2, 2, 2.40377E-01); 
%   3 <- g'
set_sigma_s(mat, m, 3, 1, 2.83340E-04);
set_sigma_s(mat, m, 3, 2, 5.24350E-02);
set_sigma_s(mat, m, 3, 3, 1.83425E-01); 
%   4 <- g'
set_sigma_s(mat, m, 4, 1, 1.46220E-06);
set_sigma_s(mat, m, 4, 2, 2.49900E-04);
set_sigma_s(mat, m, 4, 3, 9.22880E-02); 
set_sigma_s(mat, m, 4, 4, 7.90769E-02);  
set_sigma_s(mat, m, 4, 5, 3.73400E-05);  
%   5 <- g'
set_sigma_s(mat, m, 5, 1, 2.06420E-08);  
set_sigma_s(mat, m, 5, 2, 1.92390E-05);  
set_sigma_s(mat, m, 5, 3, 6.93650E-03);  
set_sigma_s(mat, m, 5, 4, 1.69990E-01);  
set_sigma_s(mat, m, 5, 5, 9.97570E-02);  
set_sigma_s(mat, m, 5, 6, 9.17420E-04);  
%   6 <- g'
set_sigma_s(mat, m, 6, 2, 2.98750E-06);  
set_sigma_s(mat, m, 6, 3, 1.07900E-03);  
set_sigma_s(mat, m, 6, 4, 2.58600E-02);  
set_sigma_s(mat, m, 6, 5, 2.06790E-01);  
set_sigma_s(mat, m, 6, 6, 3.16774E-01);   
set_sigma_s(mat, m, 6, 7, 4.97930E-02);  
%   7 <- g'
set_sigma_s(mat, m, 7, 2, 4.21400E-07);  
set_sigma_s(mat, m, 7, 3, 2.05430E-04);  
set_sigma_s(mat, m, 7, 4, 4.92560E-03);  
set_sigma_s(mat, m, 7, 5, 2.44780E-02);  
set_sigma_s(mat, m, 7, 6, 2.38760E-01);   
set_sigma_s(mat, m, 7, 7, 1.09910E+00);  

% ---------------------------
% (f) guide tube
% ---------------------------
m = 6;
% Transport cross section
set_sigma_t(mat, m, 1, 1.26032E-01);
set_sigma_t(mat, m, 2, 2.93160E-01);
set_sigma_t(mat, m, 3, 2.84240E-01);
set_sigma_t(mat, m, 4, 2.80960E-01);
set_sigma_t(mat, m, 5, 3.34440E-01);
set_sigma_t(mat, m, 6, 5.65640E-01);
set_sigma_t(mat, m, 7, 1.17215E+00);
% Absorption cross section
% set_sigma_a(mat, m, 1, 5.11320E-04);
% set_sigma_a(mat, m, 2, 7.58010E-05);
% set_sigma_a(mat, m, 3, 3.15720E-04);
% set_sigma_a(mat, m, 4, 1.15820E-03);
% set_sigma_a(mat, m, 5, 3.39750E-03);
% set_sigma_a(mat, m, 6, 9.18780E-03);
% set_sigma_a(mat, m, 7, 2.32420E-02);
% Scattering
%   1 <- g'
set_sigma_s(mat, m, 1, 1, 6.61659E-02);
%   2 <- g'
set_sigma_s(mat, m, 2, 1, 5.90700E-02);
set_sigma_s(mat, m, 2, 2, 2.40377E-01); 
%   3 <- g'
set_sigma_s(mat, m, 3, 1, 2.83340E-04);
set_sigma_s(mat, m, 3, 2, 5.24350E-02);
set_sigma_s(mat, m, 3, 3, 1.83297E-01); 
%   4 <- g'
set_sigma_s(mat, m, 4, 1, 1.46220E-06);
set_sigma_s(mat, m, 4, 2, 2.49900E-04);
set_sigma_s(mat, m, 4, 3, 9.23970E-02); 
set_sigma_s(mat, m, 4, 4, 7.88511E-02);  
set_sigma_s(mat, m, 4, 5, 3.73330E-05);  
%   5 <- g'
set_sigma_s(mat, m, 5, 1, 2.06420E-08);  
set_sigma_s(mat, m, 5, 2, 1.92390E-05);  
set_sigma_s(mat, m, 5, 3, 6.94460E-03);  
set_sigma_s(mat, m, 5, 4, 1.70140E-01);  
set_sigma_s(mat, m, 5, 5, 9.97372E-02);  
set_sigma_s(mat, m, 5, 6, 9.17260E-04);  
%   6 <- g'
set_sigma_s(mat, m, 6, 2, 2.98750E-06);  
set_sigma_s(mat, m, 6, 3, 1.07900E-03);  
set_sigma_s(mat, m, 6, 4, 2.58600E-02);  
set_sigma_s(mat, m, 6, 5, 2.06790E-01);  
set_sigma_s(mat, m, 6, 6, 3.16774E-01);   
set_sigma_s(mat, m, 6, 7, 4.97930E-02);  
%   7 <- g'
set_sigma_s(mat, m, 7, 2, 4.21400E-07);  
set_sigma_s(mat, m, 7, 3, 2.05430E-04);  
set_sigma_s(mat, m, 7, 4, 4.92560E-03);  
set_sigma_s(mat, m, 7, 5, 2.44780E-02);  
set_sigma_s(mat, m, 7, 6, 2.38760E-01);   
set_sigma_s(mat, m, 7, 7, 1.09910E+00); 

% ---------------------------
% (g) moderator
% ---------------------------
m = 7;
% Transport cross section
set_sigma_t(mat, m, 1, 1.59206E-01);
set_sigma_t(mat, m, 2, 4.12970E-01);
set_sigma_t(mat, m, 3, 5.90310E-01);
set_sigma_t(mat, m, 4, 5.84350E-01);
set_sigma_t(mat, m, 5, 7.18000E-01);
set_sigma_t(mat, m, 6, 1.25445E+00);
set_sigma_t(mat, m, 7, 2.65038E+00);
% Absorption cross section
% set_sigma_a(mat, m, 1, 6.01050E-04);
% set_sigma_a(mat, m, 2, 1.57930E-05);
% set_sigma_a(mat, m, 3, 3.37160E-04);
% set_sigma_a(mat, m, 4, 1.94060E-03);
% set_sigma_a(mat, m, 5, 5.74160E-03);
% set_sigma_a(mat, m, 6, 1.50010E-02);
% set_sigma_a(mat, m, 7, 3.72390E-02);
% Scattering
%   1 <- g'
set_sigma_s(mat, m, 1, 1, 4.44777E-02);
%   2 <- g
set_sigma_s(mat, m, 2, 1, 1.13400E-01);
set_sigma_s(mat, m, 2, 2, 2.82334E-01); 
%   3 <- g'
set_sigma_s(mat, m, 3, 1, 7.23470E-04);
set_sigma_s(mat, m, 3, 2, 1.29940E-01);
set_sigma_s(mat, m, 3, 3, 3.45256E-01); 
%   4 <- g'
set_sigma_s(mat, m, 4, 1, 3.74990E-06);
set_sigma_s(mat, m, 4, 2, 6.23400E-04);
set_sigma_s(mat, m, 4, 3, 2.24570E-01); 
set_sigma_s(mat, m, 4, 4, 9.10284E-02);  
set_sigma_s(mat, m, 4, 5, 7.14370E-05);  
%   5 <- g'
set_sigma_s(mat, m, 5, 1, 5.31840E-08);  
set_sigma_s(mat, m, 5, 2, 4.80020E-05);  
set_sigma_s(mat, m, 5, 3, 1.69990E-02);  
set_sigma_s(mat, m, 5, 4, 4.15510E-01);  
set_sigma_s(mat, m, 5, 5, 1.39138E-01);  
set_sigma_s(mat, m, 5, 6, 2.21570E-03);  
%   6 <- g'
set_sigma_s(mat, m, 6, 2, 7.44860E-06);  
set_sigma_s(mat, m, 6, 3, 2.64430E-03);  
set_sigma_s(mat, m, 6, 4, 6.37320E-02);  
set_sigma_s(mat, m, 6, 5, 5.11820E-01);  
set_sigma_s(mat, m, 6, 6, 6.99913E-01);   
set_sigma_s(mat, m, 6, 7, 1.32440E-01);  
%   7 <- g'
set_sigma_s(mat, m, 7, 2, 1.04550E-06);  
set_sigma_s(mat, m, 7, 3, 5.03440E-04);  
set_sigma_s(mat, m, 7, 4, 1.21390E-02);  
set_sigma_s(mat, m, 7, 5, 6.12290E-02);  
set_sigma_s(mat, m, 7, 6, 5.37320E-01);   
set_sigma_s(mat, m, 7, 7, 2.48070E+00);   

end