%% Load files
load('build/results/analytical.txt');
load('build/results/incremental0.txt');
load('build/results/incrementalResidual0.txt');
load('build/results/incremental1.txt');
load('build/results/incrementalResidual1.txt');
load('build/results/incremental2.txt');
load('build/results/incrementalResidual2.txt');
load('build/results/incrementalAbove.txt');

%% Parameters
E = 70000000000;
A = 0.01;
a = 2;
b = 1;
l0 = sqrt(a^2 + b^2);

qcr = sqrt(3)/9 * A*E*b^3/l0^3;

p = linspace(-0.5,0.5, 5);
q = E*A*b^3 / (l0^3*qcr) * p;

% Display results
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 25]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(analytical(1,:), analytical(2,:), 'k', 'linewidth', 1.5);
plot(incremental0(1,:), incremental0(2,:), '-ro', 'LineWidth',2);
plot(incremental1(1,:), incremental1(2,:), '-mo', 'LineWidth',2);
plot(incremental2(1,:), incremental2(2,:), '-o', 'LineWidth',2);

leg = legend('Exact solution','5 steps', '10 steps', '20 steps',...
    'Location','northeast');
set(leg,'Interpreter','latex')



xlabel('$p$ [m]','Interpreter','latex','FontSize',28);
ylabel('$\lambda$ [-]','Interpreter','latex','FontSize',28);
xlim([0,0.5]);
ylim([0,1.2]);

%% Display residuals
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 15]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

error0 = abs(incrementalResidual0./qcr);
error1 = abs(incrementalResidual1./qcr);
error2 = abs(incrementalResidual2./qcr);

plot(incremental0(2,:), error0, '-ro', 'LineWidth',2);
plot(incremental1(2,:), error1, '-mo', 'LineWidth',2);
plot(incremental2(2,:), error2, '-o', 'LineWidth',2);

leg = legend('5 steps', '10 steps', '20 steps',...
    'Location','northeast');
set(leg,'Interpreter','latex')



xlabel('$\lambda$ [-]','Interpreter','latex','FontSize',28);
ylabel('$OOB/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
%xlim([0,0.5]);
%ylim([0,1.2]);

%% Above qcr
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 15]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(analytical(1,:), analytical(2,:)./1.5, 'k', 'linewidth', 1.5);
plot(incrementalAbove(1,:), incrementalAbove(2,:), '-ro', 'LineWidth',2);

leg = legend('Exact solution','Incremental solution',...
    'Location','northeast');
set(leg,'Interpreter','latex')


xlabel('$p$ [m]','Interpreter','latex','FontSize',28);
ylabel('$\lambda$ [-]','Interpreter','latex','FontSize',28);
xlim([-1,2]);
ylim([-1,1.5]);

%% Finite difference error (old version of the code)
%{
h = logspace(-4, -1, 10);

p = 0.3;

beta = E*A/(l0^3);
exact = beta * (3*p^2 - 6*b*p + 2*b^2);
FD = (OOBmatlab(beta, b, p+h/2) - OOBmatlab(beta, b, p-h/2))./h;

error = (FD - exact) ./ exact;

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 15]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(h, error, '-o', 'linewidth', 2);

xlabel('$h$ [m]','Interpreter','latex','FontSize',28);
ylabel('Relative error [-]','Interpreter','latex','FontSize',28);
set(gca,'XScale', 'log', 'YScale', 'log')

%xlim([-1,2]);
%ylim([-1,1.5]);
%}






