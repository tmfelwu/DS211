function lambda = intersection (dc, d, delta)
    a = d'*d;
    b = 2*dc'* d;
    c = dc' * dc - delta ^ 2;
    lambda = ( -b + ( b^2 - 4 * a*c)^0.5)/ (2*a);
end