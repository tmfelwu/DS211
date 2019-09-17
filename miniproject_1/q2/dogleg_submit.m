clc;clear;
tic;
f = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;


grad_f = @(x)  [-400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
                         200*(x(2)-x(1)^2) ];
                     

hessian_f = @(x) [2 - 400*x(2) + 1200*x(1)^2,  -400*x(1);
                    -400*x(1), 200   ];
                
              
x_k = [-1.5;1.5]; %initial starting
del_max = 8;           
delta = 2;              
neta = 1/8;              
max_iter = 300;             
iter = 0;               
tol = 1E-10; %stopping tolerance           
change = inf;           

x_all = x_k; 

while change > tol && iter <= max_iter
   
    B = hessian_f(x_k);
    g = grad_f(x_k);
    f_k = f(x_k);
    p_B = -B\g; 
    
    % Model equation
    m = @(p) f_k + g'*p + 1/2*p'*B*p;
    
    if norm(p_B) <= delta
        p_k = p_B;
    else
        p_U = -g'*g/(g'*B*g)*g;
        if norm(p_U) >= delta
            p_k = delta*p_U/np_U;
        else
            % Find tau as positive root of polynomial
            p_C = p_B - p_U;
            coeffs = [norm(p_C)^2, 2*p_C'*p_U, (norm(p_U))^2-delta^2];
            tau = max(roots(coeffs));
            p_k = p_U + tau*p_C;
        end
    end
    
    rho = (f_k - f(x_k + p_k))/(f_k - m(p_k)); %reduction ratio
    if rho < 1/4 
        delta = 1/4*delta;
    else 
        if rho > 3/4 && norm(p_k) == delta
            delta = min(2*delta,del_max);
        end
    end
    if rho > neta
        x_new = x_k + p_k;
        change = abs(f(x_new)-f_k);
    else
        x_new = x_k;
        change = inf;
    end
    
    x_k = x_new;
    x_all(:,end+1) = x_k;
    iter = iter + 1;
end

if iter <= max_iter
    disp(['Dogleg method converged in ' int2str(iter) ' iterations.'])
else
    disp('Dogleg method failed to converge in the maximum amount of iterations')
end

%Plot path taken by Dogleg method

for i = 2:length(x_all)
    plot(x_all(1,i-1:i),x_all(2,i-1:i),'linewidth',3)
    hold on
    plot(x_all(1,i-1),x_all(2,i-1),'r*','linewidth',3)
end
    
   xlabel('X')
   ylabel('Y')
   title('Convergence plot of Rosenbrock Function by Trust Region Dogleg Method')



toc;
