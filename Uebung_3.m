% Aufgabe 4
f_print = @(x0, x_opt, fnct) fprintf("Minimum from [%0.2f, %0.2f]:\t %0.4f @ [%0.2f, %0.2f]\n", x0(1), x0(2), fnct(x_opt), x_opt(1), x_opt(2));
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
f_print(x0, x_opt, f);
fprintf("-------------------------------------------\n");
fprintf("-------------------------------------------\n");

% d)
for xi = x
    for yi = y
        x0 = [xi, yi];
        x_opt = fminsearch(f, x0);
        %f_print(x0, x_opt, f);
    end
end

fprintf("-------------------------------------------\n");
fprintf("-------------------------------------------\n");

% Aufgabe 5
f_rosen_orig = @(x, y) 100*(y - x.^2).^2 + (1 - x).^2;
f_rosen = @(x) f_rosen_orig(x(1), x(2));

% a)
for xi = x
    for yi = y
        x0 = [xi, yi];
        x_opt = fminsearch(f_rosen, x0);
        %f_print(x0, x_opt, f_rosen);
    end
end

fprintf("-------------------------------------------\n");
fprintf("-------------------------------------------\n");

% b)
x0 = [0.25, -0.5];
options = optimset('Display','iter','PlotFcns',@optimplotfval);
fminsearch(f_rosen, x0, options);

%c
%x0 = zeros(50,1);
x0 = 0.8 * ones(50,1);
dims = numel(x0);
f_rosen_multi_dim = @rosen_multi_dim;
x_opt = fminsearch(f_rosen_multi_dim, x0);

fprintf( ['Solution of %d-dimensional Rosenbrock: [', repmat('%0.4f ', 1, dims), ']\n'], dims, x_opt);

function f_val = rosen_multi_dim(x)
    f_val = 0;
    dims = numel(x);
    for i = 1:(dims-1)
        f_val = f_val + ( 100*(x(i+1) - x(i).^2).^2 + (1 - x(i)).^2 );
    end
end