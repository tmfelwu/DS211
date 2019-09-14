%% BFGS 
clear all;
clc;
syms x1 x2;
f = symfun( 100*(x2 - x1^2)^2 + ( 1- x1)^2 , [x1 x2]);
global iterates;
norms = [];
BFGS(f, [1,-0.5]', eye(2,2), 0.1)

%% Plot BFGS solution path
hold on
x1 = -2:0.01:1.5;
x2 = -1:0.01:2;
[X1,X2] = meshgrid(x1,x2);
z = 100 * ( X2 - X1.^2).^2 + (1 - X1).^2;
contour(X1,X2,z, 'LevelList', [-50:2:50 , 50:10:200])
plot(iterates(1,:), iterates(2,:),'-*')
hold off