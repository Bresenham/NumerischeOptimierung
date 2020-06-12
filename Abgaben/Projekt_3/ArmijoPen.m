% Fach: Numerische Optimierung
% Projekt 3 - Penalty-Verfahren & SQP-Verfahren
%
% Autor: Maximilian Gaul
% Date: 12.06.2020
%---------------------------------------------------------------

function sig = ArmijoPen(f,r,x,d,si_start)
% Schrittweitensteuerung mit der Armijo-Regel mit Aufweitung
% Input-Parameter:
% f: Function handle fuer Funktion und Gradient
% r: Penalty-Parameter
% x: Aktuelle Naeherung
% d: Aktuelle Suchrichtung
% si_start: Startschrittweite (meist letzte Schrittweite oder 1)

sig=si_start;beta=0.9;c=0.01;
[f0,gf0]=f(x,r);
[f1,gf1]=f(x+sig*d,r);

if(f1<f0+c*sig*gf0*d) % Weitere Aufweitung noetig: vergroessere sigma
    while (f1<f0+0.01*sig*gf0*d)
        sig=sig/beta;
        [f0,gf0]=f(x,r);
        [f1,gf1]=f(x+sig*d,r);
    end
    sig=sig*beta; % Korrektur des letzten Schrittes
else     % Armijo-Bedingung nicht erfuellt: verkleinere sigma
    while (f1>f0+c*sig*gf0*d)
        sig=sig*beta;
        [f0,gf0]=f(x,r);
        [f1,gf1]=f(x+sig*d,r);
    end
end
