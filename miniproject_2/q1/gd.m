function [t, y] = gd
% Generate Data
    t = 0:0.01:5;
    mean = 0;
    std = 2;
    r = 2.7;
    y = exp(-r*t) + normrnd( mean , std,size(t))  ;
end