%% Parameters
E = 70000000000;
A = 0.01;
a = 2;
b = 1;
l0 = sqrt(a^2 + b^2);

qcr = sqrt(3)/9 * A*E*b^3/l0^3;

p = linspace(-0.5,0.5, 5);
q = E*A*b^3 / (l0^3*qcr) * p;

%% To check some results
E = 70000000000;
A = 0.01;
a = 2;
b = 1;
l0 = sqrt(a^2 + b^2);

qef = 1.2*2*12049281;
lambda = 1.51856;

u = 0.50505;
v = 0.594243;
P = 0;
Q = lambda*qef;

l = l0 + u - v;

E*A/(2*l0^3) * (2*u*(u-b)*(u-2*b) + (l-l0)*(l+l0)*l) - P

- E*A/(2*l0^3) * (l-l0)*(l+l0)*l - Q

%% For linear approximation

p = linspace(-0.5,2.5, 100);

qLin = 2*E*A*b^3 / (l0^3) .* p;
qExact = 2*E*A/(2*l0^3) * p.*(p-b).*(p-2*b);

sinAlpha = b/l0;
qEps = 2*E*A* (p/l0 - sinAlpha) .* (sqrt(1 - 2*p/l0*sinAlpha + p.^2/l0^2) - 1);


%error = (qLin - qExact)./qExact;

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 45 20]);
set(gca, 'fontsize',28);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(p, qExact, 'linewidth', 2);
plot(p, qEps, 'm','linewidth', 2);
plot(p(10:25), qLin(10:25), 'r', 'linewidth', 2);


xlabel('$p$ [m]','Interpreter','latex','FontSize',28);
ylabel('Force [N]','Interpreter','latex','FontSize',28);

leg = legend('$\mathbf{S}=E\mathbf{E}^{GL}$','$\mathbf{\sigma} = E\mathbf{\varepsilon}$', 'Linear solution',...
    'Location','southeast');
set(leg,'Interpreter','latex')







