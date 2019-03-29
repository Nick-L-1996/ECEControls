%Nicholas Lanotte Lab 1

%Parameters
g = 9.8;
L = 11.2;
Mp = 381;
Ip = 0.00616;
Mw = 36;
Icmw = 0.00000746;
rw = 2.1;
R = 4.4;
kb = 0.444;
kt = 0.470;

q1 = Icmw +rw*(Mp*(L+rw)+Mw*rw);
q2 = Icmw*(Ip+(L^2)*Mp)+(rw^2)*(Ip*(Mp+Mw)+(L^2)*Mp*Mw);
q3 = Ip+L*Mp*(L+rw);

%State Space Matrices

A = [0 1 0 0;(g*L*Mp((Mp+Mw)*(rw^2)+Icmw))/q2 -q1*kb*kt/(q2*R) 0 -q1*kb*kt/(q2*R*rw);
    0 0 0 1; g*(L^2)*(Mp^2)*(rw^2)/q2 -q3*kb*kt*rw/(q2*R) 0 -q3*kb*kt/(q2*R);];

B = [0;-q1*kt/(q2*R);0;-q3*kt*rw/(q2*R);];

C = [1 0 1 0];

D = 0;

%Question 1

    