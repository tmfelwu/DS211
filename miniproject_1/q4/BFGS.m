function [ iterates, gradients ] = BFGS( f , x0, H0, epsi)
    %BFGS Broydon Fletcher Goldfarb Shanno
    %   @param {f} Symbolic Objective function
    %   @param {x0} Starting point
    %   @param {H0} initial hessian approximation
    %   @param {epsi} Tolerance for the norm of the gradient
    k = 0;
    H = cell(2,1); x = H;
    n = length(symvar(f));
    H{1} = H0;
    x{1} = x0;
    iterates = [];
    gradients = [];
    while norm(gradf(f,x{1}')) > epsi  
    %while k < 50
        iterates = [ iterates , x{1}];
        gradients = [gradients, norm(gradf(f,x{1}'))];
        
        % Computing the search direction 
        p = - H{1} * gradf(f, x{1}');
        
        % Backtracking
        alpha = 1;
        c = 1/3;
        
        possStep = x{1} + alpha * p;
        while f(possStep(1), possStep(2)) >= f(x{1}(1),x{1}(2)) + c * alpha * gradf(f,x{1}')' * p  
            alpha = alpha/3;
            possStep = x{1} + alpha * p;
        end 
        
        x{2} = x{1} + alpha * p;

        s = x{2} - x{1};

        y = gradf(f,x{2}') - gradf(f,x{1}');
        
        rho = 1/(y'*s);
        H{2} = ( eye(n) - rho * s * y')*H{1}*(eye(n) - rho * y * s') + rho*s*s';

        H{1} = H{2};
        x{1} = x{2};
        k = k +1
    end 
end