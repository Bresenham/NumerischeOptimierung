%
% Autor: Maximilian Gaul
% Date: 01.05.2020
%---------------------------------------------------------

% Aufgabe 1
h = @(x) exp(-x) + 0.5*x.^2;
start_intrvl = [0,1];
min = Bisektion(h, start_intrvl);
fprintf("Minimum %0.4f mit Funktionswert %0.4f\n", min, h(min));