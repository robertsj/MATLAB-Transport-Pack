%> @file  flare_mat.m
%> @brief Three material library for FLARE study.
%>
%> The materials defined below correspond to 0, 15, and 30 MWd/kg burnups 
%> of a WH 17x17 assembly with 92 IFBA pins at 4.25% enrichment.  The
%> bundle was specified as part of a 22.251 assignment.
%> 
%> All materials are set with 
%>   set_xyz(object, material_id, group_to, group_from)
%> with the last argument applying to scattering only.
%>
function mat = flare_mat(which)

if which == 1
    
    mat = Materials(2, 2);
    % fresh
    set_diff_coef(mat, 1, 1,    1.4493e+00);
    set_diff_coef(mat, 1, 2,    3.8070e-01);
    set_sigma_t(mat, 1, 1,      0.2300);
    set_sigma_t(mat, 1, 2,      0.8756);
    set_sigma_s(mat, 1, 1, 1,   0.2300 - 2.5000e-02);
    set_sigma_s(mat, 1, 2, 1,   1.5100e-02);
    set_sigma_s(mat, 1, 2, 2,   0.8756 - 1.0420e-01);
    set_nu_sigma_f(mat, 1, 1,   7.9000e-03);
    set_nu_sigma_f(mat, 1, 2,   1.6920e-01);
    set_chi(mat, 1, 1,          1.0);

    % once
    set_diff_coef(mat, 2, 1,    1.4479e+00);
    set_diff_coef(mat, 2, 2,    3.7080e-01);
    set_sigma_t(mat, 2, 1,      0.2302);
    set_sigma_t(mat, 2, 2,      0.8990);
    set_sigma_s(mat, 2, 1, 1,   0.2302 - 2.5800e-02);
    set_sigma_s(mat, 2, 2, 1,   1.4800e-02);
    set_sigma_s(mat, 2, 2, 2,   0.8990 - 1.2000e-01);
    set_nu_sigma_f(mat, 2, 1,   7.9000e-03);
    set_nu_sigma_f(mat, 2, 2,   1.6920e-01);
    set_chi(mat, 2, 1,          1.0);

    % twice
    set_diff_coef(mat, 3, 1,    1.4479e+00);
    set_diff_coef(mat, 3, 2,    3.7080e-01);
    set_sigma_t(mat, 3, 1,      0.2300);
    set_sigma_t(mat, 3, 2,      0.9068);
    set_sigma_s(mat, 3, 1, 1,   0.2300 - 2.6200e-02);
    set_sigma_s(mat, 3, 2, 1,   1.4700e-02);
    set_sigma_s(mat, 3, 2, 2,   0.9068 - 1.1910e-01);
    set_nu_sigma_f(mat, 3, 1,   6.0000e-03);
    set_nu_sigma_f(mat, 3, 2,   1.6250e-01);
    set_chi(mat, 3, 1,          1.0);

else
  
    mat = Materials(8, 2);    
    D1 = [1.4360000  1.4366000  0.0  1.4389000  1.4381000  1.4385000  1.4389000  1.4393000];
    D2 = [0.3635000  0.3636000  0.0  0.3638000  0.3665000  0.3665000  0.3679000  0.3680000];
    R1 = [0.0272582  0.0272995  0.0  0.0274640  0.0272930  0.0273240  0.0272900  0.0273210];
    A2 = [0.0750580  0.0784360  0.0  0.0914080  0.0848280  0.0873140  0.0880240  0.0905100];
    F1 = [0.0058708  0.0061908  0.0  0.0074527  0.0061908  0.0064285  0.0061908  0.0064285];
    F2 = [0.0960670  0.1035800  0.0  0.1323600  0.1035800  0.1091100  0.1035800  0.1091100];
    S12 =[0.0177540  0.0176210  0.0  0.0171010  0.0172900  0.0171920  0.0171250  0.0170270]; 

    for m = 1:8
        set_diff_coef(mat, m, 1,    D1(m)               );
        set_diff_coef(mat, m, 2,    D2(m)               );
        set_sigma_t(mat, m, 1,      1/3/D1(m)           );
        set_sigma_t(mat, m, 2,      1/3/D2(m)           );
        set_sigma_s(mat, m, 1, 1,   1/3/D1(m) - R1(m)	);
        set_sigma_s(mat, m, 2, 1,   S12(m)             	);
        set_sigma_s(mat, m, 2, 2,   1/3/D2(m) - A2(m)   );
        set_nu_sigma_f(mat, m, 1,   F1(m)               );
        set_nu_sigma_f(mat, m, 2,   F2(m)               );
        set_chi(mat, m, 1,          1.0                 );
    end

end