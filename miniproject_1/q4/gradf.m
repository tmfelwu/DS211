function [ grad ] = gradf( f , x )
    %Grdient of f at x
    grad = double(subs(gradient(f), symvar(f), x ));
end