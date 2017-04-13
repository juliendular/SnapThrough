%% Load files
load('build/results/analytical.txt');
load('build/results/AL.txt');
load('build/results/info1.txt');

%phi0 = info1;
%phi1 = info1;
%phi2 = info1;
%phi3 = info1;

%%
% Parameters
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
%axis equal
plot(analytical(1,:), analytical(2,:)/2, 'k', 'linewidth', 1.5);
plot(AL(1,:), AL(2,:), '-ro', 'LineWidth',2);
scatter(info1(1,:), info1(2,:), 'b', 'linewidth', 3);
%scatter(phi0(1,:), phi0(2,:), 'b', 'linewidth', 3);
%scatter(phi1(1,:), phi1(2,:), 'm', 'linewidth', 3);
%scatter(phi2(1,:), phi2(2,:), 'r', 'linewidth', 3);
%scatter(phi3(1,:), phi3(2,:), 'o', 'linewidth', 3);

leg = legend('Exact relation','$\psi=0$', '$\psi=0.5/q_{ef}$', '$\psi=1/q_{ef}$', '$\psi=2/q_{ef}$',...
    'Location','southeast');
set(leg,'Interpreter','latex')



xlabel('$p$ [m]','Interpreter','latex','FontSize',28);
ylabel('$\lambda$ [-]','Interpreter','latex','FontSize',28);
%xlim([-0.1,0.6]);
%ylim([0,0.6]);
