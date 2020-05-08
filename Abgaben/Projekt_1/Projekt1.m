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
% Startintervall
h0 = [0,1];
min = Bisektion(h, h0);
fprintf("Minimum x = %0.8f mit f(x) = %0.4f\n", min, h(min));

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
[min, k] = Mutation(f, f0);
fprintf("Minimum x = [ %s] mit f(x) = %0.8f after %d steps.\n", sprintf("%0.4f ", min), f(min), k);

fprintf("Funktion g:\n");
[min, k] = Mutation(g, g0);
fprintf("Minimum x = [ %s] mit g(x) = %0.8f after %d steps.\n", sprintf("%0.4f ", min), g(min), k);

% Aufgabe 4
% Siehe Erläuterung in Ausarbeitung

fprintf("Aufgabe 4\n");

% Berechnen der Durchschnittlichen Laufzeit über 'count' Durchläufe
% der verschiedenen Algorithmen.

% Mutation und fminsearch nehmen kein Intervall als Startwert
% (so wie bei Bisektion), sondern ein Skalar. Daher hier die zusätzliche
% Variable 'h_startwert'.
count = 200;
sums = zeros(count, 1);
options.TolX = 1e-6;
h_startwert = 0.8;
for i=1:count
    tic;
    %min = fminsearch(h, h_startwert, options);
    %min = Bisektion(h, h0);
    min = Mutation(h, h_startwert);
    elapsed = toc;
    sums(i) = elapsed;
end
fprintf("Durchschnittliche Laufzeit über %d Durchläufe: %0.6fms\n", count, sum(sums) * 1000 / count);

% Aufgabe 5
% Siehe Datei 'fminsearch2.m' und Erläuterung in Ausarbeitung

fprintf("Aufgabe 5\n");

% Berechnen der Durchschnittlichen Laufzeit über 'count' Durchläufe
options.TolX = 1e-6;
for i=1:count
    tic;
    min = fminsearch2(f, f0, options);
    elapsed = toc;
    sums(i) = elapsed;
end
fprintf("Durchschnittliche Laufzeit über %d Durchläufe: %0.6fms\n", count, sum(sums) * 1000 / count);

% Aufgabe 6
% Siehe Erläuterung in Ausarbeitung

% Aufgabe 7
% Siehe Erläuterung in Ausarbeitung
