%% Load files
load('build/results/analytical.txt');
load('build/results/AL.txt');


%% Parameters
E = 70000000000;
A = 0.01;
a = 2;
b = 1;
l0 = sqrt(a^2 + b^2);

qcr = sqrt(3)/9 * A*E*b^3/l0^3;

p = linspace(-0.5,0.5, 5);
q = E*A*b^3 / (l0^3*qcr) * p;

%% Display results
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 25]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on
plot(analytical(1,:), analytical(2,:)/2, 'k', 'linewidth', 1.5);
plot(AL(1,:), AL(2,:), '-bo', 'LineWidth',2);

leg = legend('Exact solution','Spherical arc-length method',...
    'Location','southeast');
set(leg,'Interpreter','latex')

xlabel('$p$ [m]','Interpreter','latex','FontSize',28);
ylabel('$\lambda$ [-]','Interpreter','latex','FontSize',28);
%xlim([-0.1,0.6]);
%ylim([-1,1.5]);




%% Three-bar truss

load('build/results/AL3.txt');


figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 25]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(AL3(1,:), AL3(3,:), 'b', 'LineWidth',2);
plot(AL3(2,:), AL3(3,:), 'm', 'LineWidth',2);

leg = legend('$p_0$','$p_1$',...
    'Location','southeast');
set(leg,'Interpreter','latex')


xlabel('$p_i$ [m]','Interpreter','latex','FontSize',28);
ylabel('$q_1/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
%ylim([-1.1,1.1]);
%xlim([0,2.5]);



%% Snap-back (or not) (from a previous version of the code, to get graphs)
%{
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 25]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on
%axis equal
%plot(analytical(1,:), analytical(2,:)/2, 'k', 'linewidth', 1.5);
plot(AL3_1(2,:), 2*AL3_1(3,:), 'b', 'LineWidth',2);
plot(AL3_2(2,:), 2*AL3_2(3,:), 'm', 'LineWidth',2);
plot(AL3_3(2,1:180), 2*AL3_3(3,1:180), 'color', [1 0.5 0], 'LineWidth',2);
plot(AL3_4(2,1:180), 2*AL3_4(3,1:180), 'color', [0.5 0.9 0], 'LineWidth',2);
plot(AL3_5(2,1:180), 2*AL3_5(3,1:180), '.', 'color', [0 0.9 0.9], 'LineWidth',2);


%max(abs(AL3(2,:)- AL3(3,:)))

leg = legend('$w=10$','$w=1$','$w=0.25$','$w=0.2$', '$w=0.17$',...
    'Location','southeast');
set(leg,'Interpreter','latex')

xlabel('$p_1$ [m]','Interpreter','latex','FontSize',28);
ylabel('$q_1/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
ylim([-1.1,1.1]);
xlim([0,2.5]);

%}







