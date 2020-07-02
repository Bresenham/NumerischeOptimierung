% Fach: Numerische Optimierung
% Projekt 4
%
% Autor: Maximilian Gaul
% Date: 02.07.2020
%---------------------------------------------------------------

% Aufgabe 1
% Siehe GlobNewton.m

% Aufgabe 2
tol = 1e-12;

x0 = [0; 0];
x1 = [-1.2; 1];

f_himmel = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
f_himmel_grad = @(x) [ 2 * (x(1).^2 + x(2) - 11) * 2 * x(1) + 2 * (x(1) + x(2).^2 - 7);
                    2 * (x(1).^2 + x(2) - 11) + 2 * (x(1) + x(2).^2 - 7) * 2 * x(2) ];
f_himmel_hessian = @(x) [   4 * ((x(1).^2) + x(2) - 11) + 8 * (x(1).^2) + 2, 4 * x(1) + 4 * x(2);
                            4 * x(1) + 4 * x(2), 4 * (x(1) + (x(2)).^2 - 7) + 8 * (x(2)).^2 ];
                
                
f_rosen = @(x) 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
f_rosen_grad = @(x) [ 400*x(1).^3 - 400*x(1)*x(2)+2*x(1)-2; 200*(x(2) - x(1).^2) ];
f_rosen_hessian = @(x) [ 800 * x(1)^2 - 400 * ( x(2)-x(1)^2 ) + 2, -400 * x(1); -400*x(1), 200 ];

ret = GlobNewton(f_himmel, f_himmel_grad, f_himmel_hessian, x0, tol);
ret = GlobNewton(f_himmel, f_himmel_grad, f_himmel_hessian, x1, tol);

ret = GlobNewton(f_rosen, f_rosen_grad, f_rosen_hessian, x0, tol);
ret = GlobNewton(f_rosen, f_rosen_grad, f_rosen_hessian, x1, tol);
