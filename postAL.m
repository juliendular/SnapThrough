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

%leg = legend('$p_0$','$p_1$',...
%    'Location','southeast');
%set(leg,'Interpreter','latex')


%xlabel('$p_i$ [m]','Interpreter','latex','FontSize',28);
%ylabel('$q_1/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
%ylim([-1.1,1.1]);
%xlim([0,2.5]);


%% Five-bar truss

load('build/results/AL5.txt');
load('build/results/AL7.txt');

%test1_10 = AL7;
%test2_10 = AL7;
%test4_10 = AL7;
%test6_10 = AL7;
%test7_10 = AL7; % 7.5/10
test10_10 = AL7;

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 25]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(AL5(1,:), 7*AL5(3,:), 'k', 'LineWidth',2);
plot(test1_10(1,:), 7*test1_10(3,:), 'b', 'LineWidth',2);
plot(test2_10(1,:), 7*test2_10(3,:), 'm', 'LineWidth',2);
plot(test4_10(1,:), 7*test4_10(3,:), 'r', 'LineWidth',2);
plot(test6_10(1,:), 7*test6_10(3,:), 'color', [1 0.5 0], 'LineWidth',2);
plot(test7_10(1,:), 7*test7_10(3,:), 'color', [0.5 0.9 0], 'LineWidth',2);
plot(test10_10(1,:), 7*test10_10(3,:), 'color', [0 0.9 0.9], 'LineWidth',2);


leg = legend('$\gamma = 0.0$','$\gamma = 0.1$','$\gamma = 0.2$', ...
    '$\gamma = 0.4$', '$\gamma = 0.6$','$\gamma = 0.75$','$\gamma = 1.0$',...
    'Location','northwest');
set(leg,'Interpreter','latex')


xlabel('$p$ [m]','Interpreter','latex','FontSize',28);
ylabel('$q/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
ylim([-1.1,6]);
xlim([0,2.5]);




%% Snap-back (or not) (from a previous version of the code, to get graphs)
load('build/results/AL3.txt');

%AL3_1 = AL3; % 2
%AL3_2 = AL3; % 0.5
%AL3_3 = AL3; % 0.25
%AL3_4 = AL3; % 0.125
%AL3_5 = AL3; % 0.08

%%{
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 25]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on
%axis equal
mult = 2/1.115;
%plot(analytical(1,:), analytical(2,:)/2, 'k', 'linewidth', 1.5);
plot(AL3_1(2,:), mult*AL3_1(3,:), 'b', 'LineWidth',2);
plot(AL3_2(2,:), mult*AL3_2(3,:), 'm', 'LineWidth',2);
plot(AL3_3(2,:), mult*AL3_3(3,:), 'color', [1 0.5 0], 'LineWidth',2);
plot(AL3_4(2,:), mult*AL3_4(3,:), 'color', [0.5 0.9 0], 'LineWidth',2);
plot(AL3_5(2,:), mult*AL3_5(3,:), 'color', [0 0.9 0.9], 'LineWidth',2);


%max(abs(AL3(2,:)- AL3(3,:)))

leg = legend('$w=2$','$w=0.5$','$w=0.25$','$w=0.125$', '$w=0.08$',...
    'Location','southeast');
set(leg,'Interpreter','latex')

xlabel('$p_1$ [m]','Interpreter','latex','FontSize',28);
ylabel('$q_1/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
ylim([-1.1,1.5]);
xlim([0,2.5]);

%}







