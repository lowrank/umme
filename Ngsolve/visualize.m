function visualize(v, f, data, str)
    h = figure(1);
    set(gcf,'Position',[100 100 500 400]);
    trisurf(f', v(1,:), v(2,:), data, 'EdgeColor', 'None');
    shading interp;colorbar;colormap jet;
    view(2);
    saveas(h, str);
end

