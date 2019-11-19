function jac = jacobian(x,t)
%JACOBIAN Summary of this function goes here
%   Detailed explanation goes here
   jac = -t.* exp(-x * t);
end

