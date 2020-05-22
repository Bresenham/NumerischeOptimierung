% Fach: Numerische Optimierung
% Projekt 2 - Quasi-Newton-Verfahren & Gauß-Newton-Verfahren
%
% Autor: Maximilian Gaul
% Date: 22.05.2020
%---------------------------------------------------------

% Aufgabe 1
% Siehe Datei 'InverseBFGS.m'

x0 = [0; -1];

% Funktion von Himmelblau zum Testen der BFGS-Methode
f_himmel = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
f_himmel_grad = @(x) [ 2 * (x(1).^2 + x(2) - 11) * 2 * x(1) + 2 * (x(1) + x(2).^2 - 7);
                    2 * (x(1).^2 + x(2) - 11) + 2 * (x(1) + x(2).^2 - 7) * 2 * x(2) ];

% Rosenbrock-Funktion zum Testen der BFGS-Methode
f_rosen = @(x) 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_rosen_grad = @(x) [ 400*x(1).^3 - 400*x(1)*x(2)+2*x(1)-2; 200*(x(2) - x(1).^2) ];

ret = InverseBFGS(f_himmel, f_himmel_grad, x0);
fprintf("InverseBFGS returned x=[%0.4f, %0.4f] with f_himmel(x)=%0.6f\n", ret(1), ret(2), f_himmel(ret));

ret = InverseBFGS(f_rosen, f_rosen_grad, x0);
fprintf("InverseBFGS returned x=[%0.4f, %0.4f] with f_rosen(x)=%0.6f\n", ret(1), ret(2), f_rosen(ret));