cmap =['#4b1a93';'#1a2293';'#1a6893';...
        '#1a9375';'#1a931e';'#85931a';...
        '#93801a';'#d9b045';'#d98d45'];
    
constraint_color = '#8B658B';
kk = figure;
kk.Position = [100 100 500 500];

hold on;
Name = 'Coef_1_1_1_1_OptResTraj';
for i = 1 : 9 
    file_name = Name + string(i) + ".mat";
    load(file_name);
    time = linspace(0, double(T), 10*10+1);
    x = polyval(x_poly_res, time);
    y = polyval(y_poly_res, time);
    z = polyval(z_poly_res, time);
    plot3(x, y, z, 'Color', cmap(i,:), 'LineWidth', 1.2);
    plot3(x(1:10:end), y(1:10:end), z(1:10:end), '.','MarkerSize', 8, 'Color', cmap(i,:), 'LineWidth', 0.2);
    plot3(x(1), y(1), z(1), 'o', 'Color', cmap(i,:),...
        'MarkerFaceColor', 'white', 'MarkerEdgeColor',cmap(i,:),'LineWidth', 1.2);
    plot3(x(end), y(end), z(end), 'pentagram', 'Color', 'r',...
        'MarkerSize', 10,'MarkerFaceColor', 'white', 'MarkerEdgeColor', '#CD2626','LineWidth', 1);
end
grid on;
box on;
axis equal;
xlim([-1.6, 0.1])
ylim([-0.8, 0.8])
zlim([-0.8, 0.8])
view([-1,-1,1.5])
%
figure();
subplot(3,1,1);
hold on;
for i = 1 : 9 
    file_name = Name + string(i) + ".mat";
    load(file_name);
    time = linspace(0, double(T), 22*5+1);
    x_vel_raw = linspace(7,0,8)'.* x_poly_res;
    x_vel = x_vel_raw(1:end-1);
    x = polyval(x_vel, time);
    ylim([0.05, 0.45])
    plot(time, x, 'Color', cmap(i,:), 'LineWidth', 1.2);
end
plot(time, vel_con_upp * ones(size(time)),'--', 'Color',constraint_color , 'LineWidth', 1.2);
plot(time, vel_con_low * ones(size(time)),'--', 'Color',constraint_color , 'LineWidth', 1.2);
box on;

subplot(3,1,2);
hold on;
for i = 1 : 9 
    file_name = Name + string(i) + ".mat";
    load(file_name);
    time = linspace(0, double(T), 22*5+1);
    y_vel_raw = linspace(7,0,8)'.* y_poly_res;
    y_vel = y_vel_raw(1:end-1);
    y = polyval(y_vel, time);
    ylim([-0.3, 0.3])
    plot(time, y , 'Color', cmap(i,:), 'LineWidth', 1.2);
end
box on;

plot(time, vel_con * ones(size(time)),'--', 'Color',constraint_color , 'LineWidth', 1.2);
plot(time, -vel_con * ones(size(time)),'--', 'Color',constraint_color , 'LineWidth', 1.2);

subplot(3,1,3);
hold on;
for i = 1 : 9 
    file_name = Name + string(i) + ".mat";
    load(file_name);
    time = linspace(0, double(T), 22*5+1);
    z_vel_raw = linspace(7,0,8)'.* z_poly_res;
    z_vel = z_vel_raw(1:end-1);
    z = polyval(z_vel, time);
    ylim([-0.3, 0.3])
    plot(time, z, 'Color', cmap(i,:), 'LineWidth', 1.2);
end

plot(time, vel_con * ones(size(time)),'--', 'Color',constraint_color , 'LineWidth', 1.2);
plot(time, -vel_con * ones(size(time)),'--', 'Color',constraint_color , 'LineWidth', 1.2);
box on;
