Name = 'NewReturnTest';
file_name = Name  + ".mat";
color = '#483D8B';
load(file_name);

x_poly_res = xy_poly_res(1:8);
y_poly_res = xy_poly_res(9:end);

time = linspace(0, double(T), 1000+1);
x = polyval(x_poly_res, time);
y = polyval(y_poly_res, time);
z = polyval(z_poly_res, time);
hold on;
points = 10;
plot3(x, y, z, 'color',color,'LineWidth', 1.2);
plot3(x(1:points:end), y(1:points:end), z(1:points:end),'.','color',color, 'MarkerSize', 10,  'LineWidth', 0.2);
plot3(x(1), y(1), z(1), 'o','MarkerEdgeColor',color,'MarkerFaceColor', 'white', 'LineWidth', 1.2);
plot3(x(end), y(end), z(end), 'pentagram', 'Color', 'r',...
        'MarkerSize', 10,'MarkerFaceColor', 'white', 'MarkerEdgeColor', '#CD2626','LineWidth', 1);
grid on;
box on;
axis equal;

xlabel('X')
ylabel('Y')
zlabel('Z')
view([-1,-1,1.5]);