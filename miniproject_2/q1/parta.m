[t,y] = gd();
fig = figure(1);
fig = plot(t,y,'*r','HandleVisibility','off')
xlabel('t')
ylabel('y')
title('Scatter plot of t and corresponding $$y=\exp (- 2.7 t)+ \mathcal{N}(0,2)$$','interpreter','latex')

% Plotting the original function without noise
hold on
x = 0:0.01:5;
y1 = exp( - 2.7 * x);
plot(x,y1,'b')
legend('$ y = e^{-2.7t}$','interpreter','latex')
exportfig(gcf, 'parta', 'Color', 'rgb')