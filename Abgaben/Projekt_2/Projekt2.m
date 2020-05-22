% Fach: Numerische Optimierung
% Projekt 2 - Quasi-Newton-Verfahren & Gauß-Newton-Verfahren
%
% Autor: Maximilian Gaul
% Date: 22.05.2020
%---------------------------------------------------------------

%--------------------Allgemeine Definitionen--------------------

% Funktion von Himmelblau
f_himmel = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
f_himmel_grad = @(x) [ 2 * (x(1).^2 + x(2) - 11) * 2 * x(1) + 2 * (x(1) + x(2).^2 - 7);
                    2 * (x(1).^2 + x(2) - 11) + 2 * (x(1) + x(2).^2 - 7) * 2 * x(2) ];

% Bazaraa-Shetty-Funktion
f_bazaraa = @(x) 100 * (x(1) - 2).^4 + (x(1) - 2 * x(2)).^2;
f_bazaraa_grad = @(x) [ 2 * ( 200 * (x(1) - 2).^3 + x(1) - 2 * x(2) ); 8 * x(2) - 4 * x(1) ];

% 2D Rosenbrock-Funktion
f_rosen = @(x) 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_rosen_grad = @(x) [ 400*x(1).^3 - 400*x(1)*x(2)+2*x(1)-2; 200*(x(2) - x(1).^2) ];

% N-dimensionale Rosenbrock-Funktion
f_rosen_mult = @(x) sum( 100 * ( x(2:length(x))-x(1:(length(x)-1)).^2 ).^2 + ( 1 - x(1:( length(x)-1) ) ).^2 );
f_rosen_mult_grad = @(x) f_rosen_mult_deriv_func(x);

%---------------------------------------------------------------

% Aufgabe 1
% Siehe Datei 'InverseBFGS.m'

% Aufgabe 2
% Siehe Erläuterung im PDF

% Aufgabe 3
% Siehe Erläuterung im PDF

% Aufgabe 4
% Siehe zusätzliche Erläuterung im PDF

fprintf("--------------------AUFGABE 4--------------------\n");

x0 = [0; -1];

ret = InverseBFGS(f_himmel, f_himmel_grad, x0);
fprintf("InverseBFGS returned x=%s with f_himmel(x)=%0.6f\n", vec2str( ret(end).x ), ret(end).f);

ret = InverseBFGS(f_bazaraa, f_bazaraa_grad, x0);
fprintf("InverseBFGS returned x=%s with f_bazaraa(x)=%0.6f\n", vec2str(ret(end).x), ret(end).f);

ret = InverseBFGS(f_rosen, f_rosen_grad, x0);
fprintf("InverseBFGS returned x=%s with f_rosen(x)=%0.6f\n", vec2str( ret(end).x ), ret(end).f);


rosenbrock_dim = 100;
x0_rosen = zeros(rosenbrock_dim, 1) - ones(rosenbrock_dim, 1);
ret = InverseBFGS(f_rosen_mult, f_rosen_mult_grad, x0_rosen);
fprintf("InverseBFGS returned x=%s with f_rosen_mult(x)=%0.6f\n", vec2str( ret(end).x ), ret(end).f);

% Setze die selben Toleranzen und Grenzen wie in 'InverseBFGS'
options = optimoptions("fminunc", "OptimalityTolerance", 1e-8, "MaxFunctionEvaluations", 1e+6, "MaxIterations", 1e+6);

ret = fminunc(f_rosen_mult, x0_rosen, options);
fprintf("fminunc returned x=%s with f_rosen_mult(x)=%0.6f\n", vec2str(ret), f_rosen_mult(ret));

ret = fminunc(f_himmel, x0, options);
fprintf("fminunc returned x=%s with f_himmel(x)=%0.6f\n", vec2str(ret), f_himmel(ret));

% Aufgabe 5
% Siehe Erläuterung im PDF

function ret = f_rosen_mult_deriv_func(x)

    start = @(x1, x2) -2 * (1 - x1) + 200 * (x2 - x1.^2) * (-2 * x1);
    middle = @(x1, x2, x3) 200 * (x2 - x1.^2) - 2 * (1 - x2) + 200 * (x3 -x2.^2) * (-2 * x2);
    ending = @(x2, x3) 200 * (x3 - x2.^2);
    
    dim = numel(x);
    
    sum = zeros(dim, 1);
    
    for i = 1:dim
        if i == 1
            sum(i) = sum(i) + start( x(i), x(i + 1) );
        elseif i == dim
            sum(i) = sum(i) + ending( x(i - 1), x(i) );
        else
            sum(i) = sum(i) + middle( x(i - 1), x(i), x(i + 1) );
        end
    end
    
    ret = sum;

end

function str = vec2str(vec)
    
    dim = numel(vec);
    ret_str = "[";
    for i = 1:dim
        if i == dim
            ret_str = append( ret_str, sprintf("%0.4f]", vec(i)) );
        else
            ret_str = append( ret_str, sprintf("%0.4f, ", vec(i)) );
        end
    end
    
    str = ret_str;
end