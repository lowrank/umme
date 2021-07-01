%% Case 1
[v, f, eps] = read_vtk('eps.vtk');

n = length(v);
[~, ~, sig] = read_vtk('sigma.vtk');
[~, ~, Jx] = read_vtk('Jx.vtk');
[~, ~, Jy] = read_vtk('Jy.vtk');

[~, ~, Er] = read_vtk('result_real.vtk');
[~, ~, Ec] = read_vtk('result_imag.vtk');

E = reshape( Er + 1j * Ec, 2, n);

[~, ~, Fr] = read_vtk('result_F_real.vtk');
[~, ~, Fc] = read_vtk('result_F_imag.vtk');

F = reshape( Fr + 1j * Fc, 2, n);

[~, ~, Gr] = read_vtk('result_G_real.vtk');
[~, ~, Gc] = read_vtk('result_G_imag.vtk');

G = reshape( Gr + 1j * Gc, 2, n);

[~, ~, Sr1] = read_vtk('Case_1_source_real.vtk');
[~, ~, Sc1] = read_vtk('Case_1_source_imag.vtk');

S1 = reshape( Sr1 + 1j * Sc1, 2, n);

[~, ~, Sr6] = read_vtk('Case_6_source_real.vtk');
[~, ~, Sc6] = read_vtk('Case_6_source_imag.vtk');

S6 = reshape( Sr6 + 1j * Sc6, 2, n);


