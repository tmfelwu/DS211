function y = model(x, t)
%MODEL Give the value of y for model parameter x and input values datax
%   Detailed explanation goes here
    y = exp ( -x .* t);
end