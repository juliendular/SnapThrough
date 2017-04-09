load('build/results/analytical.txt');
load('build/results/incremental.txt');
load('build/results/NR.txt');
load('build/results/AL.txt');


%%

figure;
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [10 10 40 20]);
set(gca, 'fontsize',24);
set(gca, 'fontname','timesnewroman');
box('on')
grid on
hold on

plot(analytical(1,:), analytical(2,:), 'linewidth', 1);
%plot(incremental(1,:), incremental(2,:), '-ro', 'LineWidth',2);
plot(NR(1,:), NR(2,:), '-mo', 'LineWidth',2);
plot(AL(1,:), AL(2,:), '-bo', 'LineWidth',2);


xlabel('$x$ [-]','Interpreter','latex','FontSize',24);
ylabel('$\lambda$ [-]','Interpreter','latex','FontSize',24);
%xlim([0,1]);
ylim([-1.1,1.1]);

%% 
lambda = linspace(-5, 5, 101);
x = linspace(-0.5, 2.5, 100);

for i=1:length(x)
   res(:,i) = 9/(2*sqrt(3)) * x(i) * (x(i)-1) * (x(i)-2) - lambda'; 
end

zero = zeros(1,length(analytical));
resIncr = 9/(2*sqrt(3)) * incremental(1,:) .* (incremental(1,:)-1) .* (incremental(1,:)-2) - incremental(2,:);

figure;
hold on
surf(x, lambda, res);
plot(analytical(1,:), analytical(2,:), 'linewidth', 1.5);
plot3(incremental(1,:), incremental(2,:), resIncr,'-ro', 'LineWidth',2);

shading flat;

