function [ f ] = evalf( f , x )
    %Grdient of f at x
    f = double(subs(f, symvar(f), x ));
end