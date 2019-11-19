load('data.mat');
x = 3:0.01:7;
out = zeros(1, length(x));
for i = 1:length(x)
    out(i) = norm(residuals(x(i),t,y));
end
[val, index] = min(out);
plot(x,out)
xlabel('r')
ylabel('$\frac{1}{2}\|res(r)\|_{2}^{2}$','Interpreter', 'latex')
hold on
plot(x(index), val, '*r')
caption = sprintf("r = %.2f ", x(index));
text(x(index), val -0.002 , caption, 'BackgroundColor', 'y')
exportfig(gcf,'MP2_P1_2', 'Color', 'rgb' )