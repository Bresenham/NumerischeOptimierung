% Fach: Numerische Optimierung
% Projekt 1 - Ableitungsfreie Methoden
%
% Autor: Maximilian Gaul
% Date: 01.05.2020
%---------------------------------------------------------

% Aufgabe 1
% Siehe Datei 'Bisektion.m'

fprintf("Aufgabe 1\n");

h = @(x) exp(-x) + 0.5*x.^2;

fprintf("Funktion h:\n");
start_intrvl = [0,1];
min = Bisektion(h, start_intrvl);
fprintf("Minimum x = %0.4f mit f(x) = %0.4f\n", min, h(min));

% Aufgabe 2
% Siehe Datei 'Mutation.m'

% Aufgabe 3
% Siehe Erläuterung in Ausarbeitung

fprintf("Aufgabe 3\n");

f = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
g = @(x) 100 * (x(1) - 2).^4 + (x(1) - 2*x(2)).^2;

f0 = [2,4];
g0 = [4,2];

fprintf("Funktion f:\n");
% min_f = Mutation(f, f0);
% fprintf("Minimum x = %0.4f mit f(x) = %0.4f\n", min_f, f(min_f));

fprintf("Funktion g:\n");
% min_g = Mutation(g, g0);
% fprintf("Minimum x = %0.4f mit g(x) = %0.4f\n", min_g, g(min_g));

% Aufgabe 4
% Siehe Erläuterung in Ausarbeitung

fprintf("Aufgabe 4\n");

tic;
options.TolX=1e-2;
options.Display='iter';
min = fminsearch2(f, f0, options);
toc;
disp(min);

% Aufgabe 5
% Siehe Datei 'fminsearch2.m' und Erläuterung in Ausarbeitung

% Aufgabe 6
% Siehe Erläuterung in Ausarbeitung

% Aufgabe 7
% Siehe Erläuterung in Ausarbeitung
