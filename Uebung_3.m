% Aufgabe 4
f_print = @(x0, x_opt) fprintf("Minimum for [%0.2f, %0.2f]:\t[%0.2f, %0.2f]\n", x0(1), x0(2), x_opt(1), x_opt(2));
f_orig = @(x,y) sin(x + y) * exp(-(x.^2 + 2*y.^2));
f = @(x) f_orig(x(1), x(2));


% a)
x = -3:0.5:3;
y = -3:0.5:3;
[X,Y] = meshgrid(x,y);
surf( X, Y, f_orig(X, Y) );

% c)
x0 = [-1,1];
x_opt = fminsearch(f, x0);
f_print(x0, x_opt);

% d)
