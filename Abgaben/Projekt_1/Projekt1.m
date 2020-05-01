% Fach: Numerische Optimierung
% Projekt 1 - Ableitungsfreie Methoden
%
% Autor: Maximilian Gaul
% Date: 01.05.2020
%---------------------------------------------------------

% Aufgabe 1
fprintf("Aufgabe 1\n");

h = @(x) exp(-x) + 0.5*x.^2;
start_intrvl = [0,1];
min = Bisektion(h, start_intrvl);

fprintf("Minimum %0.4f mit Funktionswert %0.4f\n", min, h(min));

% Aufgabe 2
fprintf("Aufgabe 2\n");

f = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
g = @(x) 100 * (x(1) - 2).^4 + (x(1) - 2*x(2)).^2;

f0 = [2,4];
g0 = [4,2];

min_f = Mutation(f, f0);
min_g = Mutation(g, g0);