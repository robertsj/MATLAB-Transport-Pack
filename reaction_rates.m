function rates = reaction_rates(mat, mesh, state, boundary)

mt = mesh_map(mesh, 'MATERIAL');

fr = zeros(number_groups(mat), 1);
ar = fr;
lr = zeros(number_groups(mat), 4);
ic = zeros(number_groups(mat), 4);

for g = 1:number_groups(mat)
    fr(g) = 0;
    ar(g) = 0;
    lr(g, :) = 0;
    phi = reshape(flux(state, g), number_cells_x(mesh), number_cells_y(mesh));
    for j = 1:number_cells_y(mesh)
        for i = 1:number_cells_x(mesh)
            fr(g) = fr(g) + ...
                phi(i, j)*mesh.dx(i)*mesh.dy(j)*mat.nu_sigma_f(mt(i, j),g);
            ar(g) = ar(g) + ...
                phi(i, j)*mesh.dx(i)*mesh.dy(j)*mat.sigma_a(mt(i, j),g);         
        end
    end
    set_group(boundary, g);
    lr(g, :) = lr(g, :) + get_leakage(boundary, Boundary.OUT);
    ic(g, :) = ic(g, :) + get_leakage(boundary, Boundary.IN);
end

rates.fr = fr;
rates.ar = ar;
rates.lr = lr;
rates.ic = ic;
rates.balance = sum(fr) / (sum(ar)+sum(sum(lr-ic)));

end