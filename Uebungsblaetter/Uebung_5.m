% Funktions-Definitionen
% f(x) = 2x1^2 - 4x1x2 + 4x2^2
f = @(x) 2*x(1).^2 - 4*x(1)*x(2) + 4*x(2).^2;
grad_f = @(x) [ 4*x(1) - 4*x(2); -4*x(1) + 8*x(2) ];

Q1 = [4, -4; -4, 8];
q1 = [0; 0];
x0 = [1; 0.5];

% f(x) = x1^2 + 20x2^2 - 4x1 - 2x2
g = @(x) x(1).^2 + 20*x(2).^2 - 4*x(1) - 2*x(2);
grad_g = @(x) [ 2*x(1) - 4; 40*x(2) - 2 ];

Q2 = [2, 0; 0, 40];
q2 = [-4; -2];

% Rosenbrock
f_rosen = @(x) 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_rosen_grad = @(x) [ 400*x(1).^3 - 400*x(1)*x(2)+2*x(1)-2; 200*(x(2) - x(1).^2) ];
f_rosen_hessian = @(x) [ 800 * x(1)^2 - 400 * ( x(2)-x(1)^2 ) + 2, -400 * x(1); -400*x(1), 200 ];

% Aufgabe 1
ret = ConjugateGradient(Q1, q1, x0);
fprintf("ConjugateGradient returned x=[%0.6f, %0.6f] with f(x) = %0.4f\n\n", ret(1), ret(2), 0.5 * ret' * Q1 * ret + q1' * ret);

ret = ConjugateGradient(Q2, q2, x0);
fprintf("ConjugateGradient returned x=[%0.6f, %0.6f] with g(x) = %0.4f\n\n", ret(1), ret(2), 0.5 * ret' * Q2 * ret + q2' * ret);

ret = pcg(Q1, -q1);
fprintf("pcg returned x=[%0.6f, %0.6f] for f(x)\n\n", ret(1), ret(2));

ret = pcg(Q2, -q2);
fprintf("pcg returned x=[%0.6f, %0.6f] for g(x)\n\n", ret(1), ret(2));

% Aufgabe 2
x0 = [1; 0.5];
ret = ConjugateGradientWolfe(f, grad_f, x0);
fprintf("ConjugateGradientWolfe returned x=[%0.6f, %0.6f] with f(x) = %0.4f\n\n", ret(1), ret(2), f(ret));

ret = ConjugateGradientWolfe(g, grad_g, x0);
fprintf("ConjugateGradientWolfe returned x=[%0.6f, %0.6f] with g(x) = %0.4f\n\n", ret(1), ret(2), g(ret));

ret = ConjugateGradientWolfe(f_rosen, f_rosen_grad, x0);
fprintf("ConjugateGradientWolfe returned x=[%0.6f, %0.6f] with f_rosen(x) = %0.4f\n\n", ret(1), ret(2), f_rosen(ret));

% Aufgabe 3
x0 = [2; 1];
ret = Newton(f_rosen, f_rosen_grad, f_rosen_hessian, x0);
fprintf("Newton returned x=[%0.6f, %0.6f] with f_rosen(x) = %0.4f\n\n", ret(1), ret(2), f_rosen(ret));

% Funktion Aufgabe 1
function ret = ConjugateGradient(Q, q, x0)
    
    k = 0;
    x = x0;
    grad_f = Q*x + q;
    d = -grad_f;

    while norm(grad_f) > 1e-8 && k < 100
        
        alpha = - ( grad_f' * d ) / ( d' * (Q*d) );
        x = x + alpha * d;
        grad_f_new = grad_f + alpha*Q*d;
        beta = ( norm(grad_f_new)^2 ) / ( norm(grad_f)^2 );
        d = -grad_f_new + beta * d;
        grad_f = grad_f_new;

        k = k + 1;
    end
    
    fprintf("ConjugateGradient took %d steps.\n", k);
    ret = x;
end

% Funktion Aufgabe 2
function ret = ConjugateGradientWolfe(f, grad, x0)
    
    k = 0;
    x = x0;
    d = -grad(x0);
    
    phi = @(a) f(x + a * d);
    phi_grad = @(a) grad(x + a * d)' * d;
    
    while norm( grad(x) ) > 1e-6 && k < 1000

        alpha = WolfePowell(phi, phi_grad);
        
        x_new = x + alpha * d;
        beta = norm( grad(x_new) )^2 / norm( grad(x) )^2;
        d = -grad(x_new) + beta * d;
        
        x = x_new;
        
        k = k + 1;
    end
    
    fprintf("ConjugateGradientWolfe took %d steps.\n", k);
    ret = x;
end

%Funktion Aufgabe 3
function ret = Newton(f, f_grad, f_hessian, x0)
    
    k = 0;
    x = x0;
    
    while norm( f_grad(x) ) > 1e-8 && k < 1000
        
        d = f_hessian(x) \ ( -f_grad(x) );
        
        x = x + d;
        
        k = k + 1;
    end
    
    fprintf("Newton took %d steps.\n", k);
    ret = x;
end

% phi = f(x + td)
% phi_grad = ∇f(x + td)^T * d
function ret = WolfePowell(phi, phi_grad)

    rho = 0.3;
    gamma = 1.5;
    alpha = 1.0;
    sigma = 1e-2;
    
    A = @(a) phi(a) <= phi(0) + sigma * a * phi_grad(0);
    W = @(a) phi_grad(a) >= rho * phi_grad(0);
    
    while(1)
        if(A(alpha) && W(alpha))
            ret = alpha;
            return
        end
        if(A(alpha) && ~W(alpha))
            alpha = gamma * alpha;
        elseif(~A(alpha))
            break
        end
    end
    
    a = 0;
    b = alpha;
    
    while(1)
        alpha = (a + b) / 2.0;
        if(A(alpha) && W(alpha))
            ret = alpha;
            return
        end
        if(A(alpha) && ~W(alpha))
            a = alpha;
        elseif(~A(alpha))
            b = alpha;
        end
    end
    
end