%% Generate the data
load('data.mat')
fun = @(x, datax ) exp( -x .* datax);
x0 = 1;
lb = [];
ub = [];
options = optimoptions('lsqcurvefit', ...
    'Display','iter-detailed', ...
    'Algorithm','levenberg-marquardt',...
    'OptimalityTolerance', 1e-20, ...
    'FunctionTolerance', 1e-20, ...
    'MaxFunctionEvaluations', 1000,...
    'MaxIterations', 200);
tic()
[x, resnorm, residual, exitFlag, output] = lsqcurvefit(fun, x0, t, y, lb, ub, options);
toc()
%% Some plotting
plot(t,y,'*r','HandleVisibility','off')
xlabel('t')
ylabel('y')
hold on
% Generating plot with solution
t_sol = 0:0.01:5;
y_sol = exp( - x * t_sol);
plot(t_sol, y_sol)
hold off