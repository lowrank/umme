function h = draw_circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot3(xunit, yunit, 10*ones(size(xunit)));
hold off