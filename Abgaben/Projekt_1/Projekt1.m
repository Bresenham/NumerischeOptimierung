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

% fprintf("Funktion h:\n");
% h0 = [0,1];
% min = Bisektion(h, h0);
% fprintf("Minimum x = %0.8f mit f(x) = %0.4f\n", min, h(min));

% Aufgabe 2
% Siehe Datei 'Mutation.m'

% Aufgabe 3
% Siehe Erläuterung in Ausarbeitung

fprintf("Aufgabe 3\n");

f = @(x) (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
g = @(x) 100 * (x(1) - 2).^4 + (x(1) - 2*x(2)).^2;

f0 = [2,4];
g0 = [4,2];

iterations = 1;
vals = zeros(iterations, 1);
for i=1:iterations
    fprintf("Funktion g:\n");
    fnct = f;
    [min, k] = Mutation(fnct, f0);
    fprintf("Minimum x = [ %s] mit g(x) = %0.8f after %d steps.\n", sprintf("%0.4f ", min), fnct(min), k);
    vals(i) = k;
end
fprintf("Average iteration count: %0.2f\n", sum(vals) / iterations);

% Aufgabe 4
% Siehe Erläuterung in Ausarbeitung

fprintf("Aufgabe 4\n");

% tic;
% options.TolX = 1e-6;
% min = fminsearch(h, 0, options);
% toc;

% tic;
% min = Bisektion(h, h0);
% toc;

% Aufgabe 5
% Siehe Datei 'fminsearch2.m' und Erläuterung in Ausarbeitung
% options.Display = 'iter';
% fminsearch2(f, f0, options);

% Aufgabe 6
% Siehe Erläuterung in Ausarbeitung

% Aufgabe 7
% Siehe Erläuterung in Ausarbeitung
