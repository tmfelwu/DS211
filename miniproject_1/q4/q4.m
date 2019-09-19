%%
clear all
func = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
grad_f = @(x)  [-400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
                         200*(x(2)-x(1)^2) ];
hessian_f = @(x) [2 - 400*x(2) + 1200*x(1)^2,  -400*x(1);
                    -400*x(1), 200   ];
global grad_f func hessian_f

%% BFGS 
clc;
syms x1 x2;
f = symfun( 100*(x2 - x1^2)^2 + ( 1- x1)^2 , [x1 x2]);
norms = [];
[iterates, gradients] = BFGS(f, [1,-0.5]', eye(2,2), 1e-8)

%% Plot BFGS solution path
hold on
x1 = -2:0.01:1.5;
x2 = -1:0.01:2;
[X1,X2] = meshgrid(x1,x2);
z = 100 * ( X2 - X1.^2).^2 + (1 - X1).^2;
xlabel('X')
ylabel('Y')
contour(X1,X2,z, 'LevelList', [-50:2:50 , 50:10:200])
plot(iterates(1,:), iterates(2,:),'-*')
hold off

%% SR1 Trust region

[iteratesSR1, gradientsSR1 ] = SR1(f, [-15,-15]', eye(2), 0.9, 1e-2, 1)
figure
hold on;
x1 = -2:0.01:1.5;
x2 = -1:0.01:2;
[X1,X2] = meshgrid(x1,x2);
z = 100 * ( X2 - X1.^2).^2 + (1 - X1).^2;
xlabel('X')
ylabel('Y')
contour(X1,X2,z, 'LevelList', [-50:2:50 , 50:10:200])
plot(iteratesSR1(1,:), iteratesSR1(2,:),'-*')
hold off;

%% Comparing Convergence part a
figure
set(gca, 'YScale', 'log')
hold on;
load NewtonGradientNorm.mat
plot(Grad_Norm, 'bd-', 'DisplayName','Newton Descent')
plot(gradients,'r*-', 'DisplayName','BFGS')
plot(gradientsSR1,'go-','DisplayName','SR1 Trust Region')
xlabel('k , Number of Iterations')
ylabel('log(||grad f||)')
hold off;
legend('show')