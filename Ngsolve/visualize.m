function visualize(filename, verbose)

if nargin < 2
    verbose = 0;
end
    [v, f, data] = read_vtk(filename);
    figure(1);
    trisurf(f', v(1,:), v(2,:), data(1:2:end), 'EdgeColor', 'None');
    shading interp;colorbar;colormap jet;
    view(2);
    if verbose
    hold on;
    draw_circle(0,0, 1);
    hold off;
    end
    figure(2);
    trisurf(f', v(1,:), v(2,:), data(2:2:end), 'EdgeColor', 'None');
    shading interp;colorbar;colormap jet;
    view(2);
    if verbose
    hold on;
    draw_circle(0,0, 1);
    hold off;
    end
end

