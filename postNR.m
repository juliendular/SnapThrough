%% Load files
load('build/results/analytical.txt');
load('build/results/NR0.txt');
load('build/results/NRResidual0.txt');
load('build/results/NR1.txt');
load('build/results/NRResidual1.txt');
load('build/results/NR2.txt');
load('build/results/NRResidual2.txt');
load('build/results/compositeLoading0.txt');
load('build/results/stresses0.txt');
load('build/results/compositeLoading1.txt');
load('build/results/stresses1.txt');

%% Parameters
E = 70000000000;
A = 0.01;
a = 2;
b = 1;
l0 = sqrt(a^2 + b^2);

qcr = 2*sqrt(3)/9 * A*E*b^3/l0^3;

p = linspace(-0.5,0.5, 5);
q = E*A*b^3 / (l0^3*qcr) * p;

% Display results
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 18]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(analytical(1,:), analytical(2,:), 'k', 'linewidth', 1.5);
plot(NR0(1,:), NR0(2,:), '-ro', 'LineWidth',2);
plot(NR1(1,:), NR1(2,:), '-mo', 'LineWidth',2);
%plot(NR2(1,:), NR2(2,:), '-bo', 'LineWidth',2);

leg = legend('Exact solution', '5 steps', '10 steps',...
    'Location','southeast');
set(leg,'Interpreter','latex')



xlabel('$p$ [m]','Interpreter','latex','FontSize',28);
ylabel('$\lambda$ [-]','Interpreter','latex','FontSize',28);
%xlim([0,0.5]);
ylim([-2,2]);

%% Iterations (from an old version of the code)
%{
load('build/results/NRiterations1.txt');
load('build/results/NRiterations2.txt');
load('build/results/NRiterations3.txt');
load('build/results/NRiterations4.txt');
load('build/results/NRiterations11.txt');

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 10]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

x = 1:5;
y(1) = 0.1;
for i=2:5
    y(i) = y(i-1)^2;
end

%plot(abs(NRiterations11(2,:)), 0.5, '-ro', 'linewidth', 1.5);


plot(0:length(NRiterations11)-1, abs(NRiterations11(2,:))/qcr, '-ro', 'linewidth', 1.5);
%plot(0:length(NRiterations2)-1, abs(NRiterations2(2,:))/qcr, '-ro', 'linewidth', 1.5);
%plot(NR0(1,:), NR0(2,:), '-ro', 'LineWidth',2);
%plot(NR1(1,:), NR1(2,:), '-mo', 'LineWidth',2);
%plot(NR2(1,:), NR2(2,:), '-ro', 'LineWidth',2);

%leg = legend('Multiplicity 1', 'Multiplicity 2',...
%    'Location','southwest');
%set(leg,'Interpreter','latex')


xlabel('Iteration [-]','Interpreter','latex','FontSize',28);
ylabel('$OOB/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
%xlim([0,20]);
%ylim([0,1.2]);

set(gca, 'YScale', 'log')
%}

%% Composite loading
time = linspace(0,4,length(compositeLoading0));

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 18]);
set(gca, 'fontsize',28);
box('on')
grid on
hold on
plot(time, compositeLoading0(2,:)/qcr, '-bo', 'LineWidth',2);
plot(time, compositeLoading1(2,:)/qcr, '-mo', 'LineWidth',2);
xlabel('Time [s]','Interpreter','latex','FontSize',28);
ylabel('$q/qcr$ [-]','Interpreter','latex','FontSize',28);
leg = legend('$q_{ef} = q_{cr}$', '$q_{ef} = 2q_{cr}$',...
    'Location','northeast');
set(leg,'Interpreter','latex')

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 18]);
set(gca, 'fontsize',28);
box('on')
grid on
hold on
plot(time, compositeLoading0(1,:), '-bo', 'LineWidth',2);
plot(time, compositeLoading1(1,:), '-mo', 'LineWidth',2);
xlabel('Time [s]','Interpreter','latex','FontSize',28);
ylabel('$p$ [m]','Interpreter','latex','FontSize',28);
leg = legend('$q_{ef} = q_{cr}$', '$q_{ef} = 2q_{cr}$',...
    'Location','northeast');
set(leg,'Interpreter','latex')

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 18]);
set(gca, 'fontsize',28);
box('on');
grid on
hold on
plot(time, stresses0/1e6, '-bo', 'LineWidth',2);
plot(time, stresses1/1e6, '-mo', 'LineWidth',2);
xlabel('Time [s]','Interpreter','latex','FontSize',28);
ylabel('$S_{11}$ [MPa]','Interpreter','latex','FontSize',28);
leg = legend('$q_{ef} = q_{cr}$', '$q_{ef} = 2q_{cr}$',...
    'Location','southeast');
set(leg,'Interpreter','latex')

%xlim([0,0.5]);
%ylim([-2,2]);





%% Display residuals
figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 15]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

error0 = abs(NRResidual0./qcr);
error1 = abs(NRResidual1./qcr);
error2 = abs(NRResidual2./qcr);

plot(NR0(2,:), error0, '-ro', 'LineWidth',2);
plot(NR1(2,:), error1, '-mo', 'LineWidth',2);
plot(NR2(2,:), error2, '-o', 'LineWidth',2);

leg = legend('5 steps', '10 steps', '20 steps',...
    'Location','northeast');
set(leg,'Interpreter','latex')


xlabel('$\lambda$ [-]','Interpreter','latex','FontSize',28);
ylabel('$OOB/q_{cr}$ [-]','Interpreter','latex','FontSize',28);
%xlim([0,0.5]);
%ylim([0,1.2]);
