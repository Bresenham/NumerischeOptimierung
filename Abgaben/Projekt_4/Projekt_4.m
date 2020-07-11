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

% Aufgabe 3
% Siehe PDF Dokument

% Aufgabe 4
% Siehe PDF Dokument

% Aufgabe 5
% Siehe PDF Dokument

% Aufgabe 6
% Siehe PDF Dokument

% Aufgabe 7
f_func = @(x) -2 * x(1) - 3 * x(2) - 4 * x(3);

f = [-2, -3, -4];
A = [1, 1, 1; 0, 3, 1; 1, 0, 0; 0, 0, 1];
b = [4; 5; 2; 3];

options = optimoptions("linprog", "OptimalityTolerance", 1e-8, "Algorithm", "dual-simplex");
ret = linprog(f, A, b, [], [], [0, 0, 0], [], options);

% Aufgabe 8
Q = [2, -2; -2, 4];
q = [-2; -6];
U = [0.5, 0.5; -1, 2; -1, 0; 0, -1];
r = [1, 2, 0, 0];
G = [];
b = [];
x0 = [-12; 13];

ret = ActiveSet(Q, q, G, U, b, r, x0);
disp(ret)

g = @(x) x(1).^2 + 2 * x(2).^2 - 2 * x(1) - 6 * x(2) - 2 * x(1) * x(2);
g1 = @(x) 0.5 * x(1) + 0.5 * x(2) - 1;
g2 = @(x) -x(1) + 2 * x(2) - 2;

x0 = [-12; 13];
conNeqG = @(x) confunNeqG(g1, g2, x);
ret = fmincon(g, x0, [], [], [], [], [0, 0], [], conNeqG);
disp(ret);

% Aufgabe 9
% Siehe PDF Dokument

% Aufgabe 10
% Siehe PDF Dokument

% Aufgabe 11
% Siehe PDF Dokument

% Aufgabe 12

f = @(x) 2 * x(1).^2 + 1.75 * x(2).^2 + 0.75 * x(3).^2 + 4500 * x(1) + 4000 * x(2) + 3500 * x(3) + 3000 * x(4) - ( 7.5 * 10^6 );
Q = [   4, 0, 0, 0;
        0, 3.5, 0, 0;
        0, 0, 1.5, 0;
        0, 0, 0, 0;];
q = [4500; 4000; 3500; 3000];
U = [   1, 0, 0, 0;
        1, 1, 0, 0;
        1, 1, 1, 0;
        -1, 0, 0, 0;
        -1, -1, 0, 0;
        -1, -1, -1, 0;
        -1, 0, 0, 0;
        0, -1, 0, 0;
        0, 0, -1, 0;
        0, 0, 0, -1];
r = [3500; 7500; 10500; -1500; -5500; -8500; 0; 0; 0; 0];
G = [1, 1, 1, 1];
b = [10000];
x0_0 = [3000; 4000; 2000; 1000];
x0_1 = [2500; 4750; 1425; 1325];
x0_2 = [1500; 6000; 2000; 500];

ret = ActiveSet(Q, q, G, U, b, r, x0_2);
disp( f(ret(end).x) );
disp(ret)

options = optimoptions("quadprog", "Algorithm", "active-set", "Display", "iter-detailed");
ret = quadprog(Q, q, U, r, G, b, [], [], x0_2, options);
disp(f(ret));
disp(ret);

% Ungleichheitsbedingungen aus Aufgabe 8 für fmincon
function [c,ceq] = confunNeqG(g1, g2, x)
    % Nonlinear inequality constraints
    c = [g1(x), g2(x)];
    % Nonlinear equality constraints
    ceq = [];
end