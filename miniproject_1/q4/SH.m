function x_star = SH( Q, b , x_cap, delta)
    global grad_f
    n = 2
    d = cell(2,1); x = d;
    x{1} = [ 0 , 0 ]';
    d{1} = -b;
    k = 1;
    while norm(grad_f(x{1})) ~= 0 & k ~= n+1
        if d{1}'* Q* d{1} <= 0 
            dc = - (b'*b)/(b'*Q*b)*b;
            lambda = intersection(dc, d{1} - dc, 2);
            x_star = x{1} + lambda * d{1};
        end
        alpha = - (d{1}'*(Q*x{1} + b))/ (d{1}'*Q* d{1});
        x{2} = x{1} + alpha * d{1};
        if norm(x{2}) > delta
            dc = - (b'*b)/(b'*Q*b)*b;
            lambda = intersection(dc, d{1} - dc, delta);
            x_star = x{1} + lambda * d{1};
        end
        beta = ((Q*x{2} + b)'*(Q*x{2}+b))/((Q*x{1} + b)'*(Q*x{1} + b));
        d{2} = - Q* x{2} - b + beta*d{1};
        k = k+1;

        x{1} = x{2};
        d{1} = d{2};
    end
    x_star = x{1};
end