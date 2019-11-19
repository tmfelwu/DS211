function jf = jacobian_f(x, t, y)
%JACOBIAN_F Summary of this function goes here
%   Detailed explanation goes here
    jf= 0.5 * sum( -2 * residuals(x,t,y) .* t .* exp(-x * t));
end