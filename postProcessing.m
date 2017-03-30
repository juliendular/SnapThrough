load('build/results/analytical.txt');
load('build/results/incremental.txt');





figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10 10 40 20]);
set(gca, 'fontsize',24);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(analytical(1,:), analytical(2,:), 'linewidth', 1);
plot(incremental(1,:), incremental(2,:), '-ro', 'LineWidth',2);

xlabel('$x$ [-]','Interpreter','latex','FontSize',24);
ylabel('$\lambda$ [-]','Interpreter','latex','FontSize',24);
xlim([0,0.5]);
ylim([0,1]);





