function sol = gn(func, jacobian_f, residuals, jacobian, t, y ,x0, max_iters )
% Gauss-Newton implementation
% @param {func} function handle to the model
% @param {jacobian_f} jacobian of the function
% @param {y}
% @param {x0} inital guess
% @param {max_iters}
    x = x0;
    i = 0;
    
    while i < max_iters
        
        res = residuals(x,t,y);
        jac = jacobian(x,t);
        Jt_J = jac * jac';
        Jt_r = jac * res';
        
        p_k = - Jt_J\Jt_r;
        % Backtracking Line Search
        alpha = 1;
        c = 1/3;
        
        possX = x + alpha * p_k;
        while 0.5*norm(residuals(possX,t,y)) >= 0.5 * norm(residuals(x,t,y)) %+ c * alpha * jacobian_f(x,t,y)' * p_k
            alpha = alpha/3;
            possX = x + alpha * p_k;
        end
        % Update the x
        x = x + alpha * p_k;
        i = i+1;
    end
    
    sol = x;
   
end