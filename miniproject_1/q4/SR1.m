function [ iterates, gradients, deltas ] = SR1( x0, B0, delta0, eps, del_max, max_iter)
    % SR1 Trust region
    % @param {x0} Starting point
    % @param {B0} Initial hessian approixmation
    % @param {delta0} trust region radius
    % @param {eps} convergance tolerance > 0
    % @param {del_max} max delta
    % @param {max_iters} max_iterations
    k = 0;
    B = cell(2,1); x = B ; delta = B;
    eta = 1e-4; % between 0 and 1e-3
    r = 0.5; % r between 0 and 1
    B{1} = B0;
    x{1} = x0;
    delta{1} = delta0;
    iterates = [];
    gradients = [];
    deltas = [];

    global func grad_f;
    while  norm(grad_f(x{1})) > eps & k < max_iter
    %while k < 10
        norm(grad_f(x{1}))
        iterates = [ iterates , x{1}];
        gradients = [ gradients, norm(grad_f(x{1}))];
        deltas = [deltas, delta{1}];
        g = grad_f(x{1});
        p_B = -inv(B{1}) * g;

        % START 1 : Taken from Q2 
        if norm(p_B) <= delta{1}
            p_k = p_B;
        else
            p_U = -g'*g/(g'*B{1}*g)*g;
            if norm(p_U) >= delta{1}
                p_k = delta{1}*p_U/norm(p_U);
            else
                % Find tau as positive root of polynomial
                p_C = p_B - p_U;
                coeffs = [norm(p_C)^2, 2*p_C'*p_U, (norm(p_U))^2-delta{1}^2];
                tau = max(roots(coeffs));
                p_k = p_U + tau*p_C;
            end
        end
        % END 1

        sk = p_k;

        y_k = grad_f( x{1} + sk) - grad_f(x{1});
        ared = func(x{1}') - func(x{1} + sk);
        pred = - ( grad_f( x{1})' * sk  + 0.5 * sk'* B{1} * sk);
        rho = ared/pred ;

        if rho > eta
            x{2} = x{1} + sk;
        else
            x{2} = x{1};
        end

        if rho > 0.75
            if norm(sk) <= 0.8 * delta{1}
                delta{2} = delta{1};
            else
                delta{2} = 2 * delta{1};
            end
        elseif rho >= 0.1 & rho <= 0.75
            delta{2} = delta{1};
        else
            delta{2} = 0.5 * delta{1};
        end


        if abs(sk'*( y_k - B{1}* sk)) >= r * norm(sk) * norm(y_k - B{1} * sk)
            temp = y_k - B{1} * sk;
            B{2} = B{1} + (temp * temp') / (temp' * sk);
        else 
            B{2} = B{1};
        end
        
        delta{1} = delta{2};
        B{1} = B{2};
        x{1} = x{2};
        k = k + 1
    end


end