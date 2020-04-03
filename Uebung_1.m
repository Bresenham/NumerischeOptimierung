% Aufgabe 2a)
x = -3:0.125:3;
y = -3:0.125:3;
[X,Y] = meshgrid(x,y);
F = (X.^2+Y.^2).*exp(sin(X+Y));
surf(X,Y,F);
figure;

hold on
% Aufgabe 2b)
x = -2:0.0125:0;
y = 1:0.0125:3;
[X,Y] = meshgrid(x,y);
F1 = 100 * (Y - X.^2).^2 + (1 - X).^2;
% F2 = -99 - 2*X + 200*Y + 201*X.^2 + 400*X*Y + 100*Y.^2;
F2 = -99-2*X+200*Y+201*X.^2+400*X.*Y+100*Y.^2;
surf(X, Y, F1);
surf(X, Y, F2);

% Aufgabe 3a)
m1 = [25 -5 15; -5 10 -15; 15 -15 29];
m2 = [2 3 0 -6; 3 7 4 8; 0 4 -2 7; -6 8 7 1];
m3 = [2 3; 3 2];
m5 = [2 4; 1 2];
disp( eig(m1) ); % positiv definit
disp( eig(m2) ); % indefinit
disp( eig(m3) );
disp( eig(m5) );

% Aufgabe 3b)
disp( eval(@sin, [pi 0.5*pi]) );

fnct = @(x) 100 * (x - x.^2).^2;
disp ( eval(fnct, [0 4]) );

function out_vals = eval(fnct, in_vals)
    out_vals = fnct(in_vals);
end