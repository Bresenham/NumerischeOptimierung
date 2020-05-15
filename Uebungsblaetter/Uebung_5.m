% Aufgabe 1
% f(x) = 2x1^2 - 4x1x2 + 4x2^2
Q1 = [4, -4; -4, 8];
q1 = [0; 0];

% f(x) = x1^2 + 20x2^2 - 4x1 - 2x2
Q2 = [0.5, 0; 0, 10];
q2 = [-4; -2];
x0 = [-5; 1];

ret = ConjugateGradient(Q2, q2, x0);
fprintf("RESULT: [%0.6f, %0.6f]\n", ret(1), ret(2));

disp( pcg(Q2, -q2) );

function ret = ConjugateGradient(Q, q, x0)
    
    k = 0;
    x = x0;
    grad_f = Q*x + q;
    d = -grad_f;

    while norm(grad_f) > 1e-14 && k < 1000
        alpha = - ( grad_f' * d ) / ( d' * (Q*d) );
        x = x + alpha * d;
        grad_f_new = grad_f + alpha*Q*d;
        beta = ( norm(grad_f_new)^2 ) / ( norm(grad_f)^2 );
        d = -grad_f_new + beta * d;
        grad_f = grad_f_new;

        k = k + 1;
    end
    
    ret = x;
end