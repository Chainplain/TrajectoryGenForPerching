Name = 'ReturnTest';
file_name = Name  + ".mat";
color = '#483D8B';
dubin_color = '#85931a';
Face_color = '#1a931e';

load(file_name);
Dubin_len = length(Dubin_path);
leng_array = 0 : Dubin_path_step : Dubin_len * Dubin_path_step-Dubin_path_step;


time = linspace(0, double(T), 10*10+1);
x = polyval(x_poly_res, time);
y = polyval(y_poly_res, time);
z = polyval(z_poly_res, time);
Dubin_z = polyval(z_poly_res, 5 * leng_array);
hold on;

draw_small_dot_ab = 10:10: length(x);

lw = 1.5;
plot3(x, y, z, 'color',color,'LineWidth', lw);
plot3(Dubin_path(:,1), Dubin_path(:,2), Dubin_z, 'color',dubin_color,'LineWidth', lw);
plot3(x(draw_small_dot_ab), y(draw_small_dot_ab), z(draw_small_dot_ab),'.','color',color, 'MarkerSize', 10,  'LineWidth', 0.2);
plot3(x(1), y(1), z(1), 'o','MarkerEdgeColor',color,'MarkerFaceColor', 'white', 'LineWidth', lw);
plot3(x(end), y(end), z(end), 'o','MarkerEdgeColor',color,'MarkerFaceColor', 'black', 'LineWidth', lw);


Name = 'ReturnTest_without_ob';
file_name = Name  + ".mat";
load(file_name);
x = polyval(x_poly_res, time);
y = polyval(y_poly_res, time);
z = polyval(z_poly_res, time);

no_ob_color = '#CD2626';
plot3(x, y, z, 'color',no_ob_color,'LineWidth', lw);
plot3(x(draw_small_dot_ab), y(draw_small_dot_ab), z(draw_small_dot_ab),'.','color',no_ob_color, 'MarkerSize', 10,  'LineWidth', 0.2);



% Define the cylinder parameters
radius = 0.3;
height = 0.6;

% Number of points around the circumference of the cylinder
numPoints = 100;

% Create the cylinder using the cylinder function
[cx, cy, cz] = cylinder(radius, numPoints);

% Scale the cylinder to the specified height
cz = cz * height;

% Plot the cylinder
surf(cx, cy, cz, 'FaceColor', Face_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2);

grid on;
box on;
axis equal;

xlabel('X')
ylabel('Y')
zlabel('Z')
view([-1,-1,1.5]);