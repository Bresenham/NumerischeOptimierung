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

x0 = [0; -1];

% ret = InverseBFGS(f_himmel, f_himmel_grad, x0);
fprintf("InverseBFGS returned x=[%0.4f, %0.4f] with f_himmel(x)=%0.6f\n", ret(1), ret(2), f_himmel(ret));

ret = InverseBFGS(f_bazaraa, f_bazaraa_grad, x0);
fprintf("InverseBFGS returned x=[%0.4f, %0.4f] with f_bazaraa(x)=%0.6f\n", ret(1), ret(2), f_bazaraa(ret));

% ret = InverseBFGS(f_rosen, f_rosen_grad, x0);
fprintf("InverseBFGS returned x=[%0.4f, %0.4f] with 2D f_rosen(x)=%0.6f\n", ret(1), ret(2), f_rosen(ret));

x0 = [0; 0; 0];
% ret = InverseBFGS(f_rosen_mult, f_rosen_mult_deriv, x0);
% fprintf("InverseBFGS returned x=[%0.4f, %0.4f, %0.4f] with f_rosen_mult(x)=%0.6f\n", ret(1), ret(2), ret(3), f_rosen_mult(ret));

function ret = f_rosen_mult_deriv_func(x)

    start = @(x1, x2) 2 * (1 - x1) * (-1) + 200 * (x2 - x1.^2) * (-2 * x1);
    middle = @(x0, x1, x2) 2 * (1 - x1) * (-1) + 200 * (x1 - x0.^2) + 200 * (x2 - x1.^2) * (-2 + x1);
    ending = @(x0, x1) 200 * (x1 - x0.^2);
    
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