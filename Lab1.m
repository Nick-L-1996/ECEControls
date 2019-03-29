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

q1 = Icmw+rw*(Mp*(L+rw)+Mw*rw);
q2 = Icmw*(Ip+(L^2)*Mp)+(rw^2)*(Ip*(Mp+Mw)+(L^2)*Mp*Mw);
q3 = Ip+L*Mp*(L+rw);

%State Space Matrices

A = [0, 1, 0, 0;
    (g*L*Mp*((Mp+Mw)*(rw^2)+Icmw))/q2,(-1)*q1*kb*kt/(q2*R),0,(-1)*q1*kb*kt/(q2*R*rw);
    0 0 0 1; 
    g*(L^2)*(Mp^2)*(rw^2)/q2,(-1)*q3*kb*kt*rw/(q2*R),0,(-1)*q3*kb*kt/(q2*R);];

B = [0;
    -q1*kt/(q2*R);
    0;
    -q3*kt*rw/(q2*R);];

C = [1 0 1 0];

D = 0;

%Question 1
StateModel = ss(A,B,C,D);
disp("A Matrix = ");
disp(A);
disp("B Matrix = ");
disp(B);
disp("C Matrix = ");
disp(C);
disp("D Matrix = ");
disp(D);

%Question 2
[a,b] = ss2tf(A,B,C,D);
disp("Transer Function");
disp(a);
disp(b);

%Question 3
TransferFunction = tf(a,b);
figure;
impulse(TransferFunction);

%Question 4
figure;
step(TransferFunction);
%need ramp response

%Question 5 
figure;
pzmap(TransferFunction);
[p,z] = pzmap(TransferFunction);
disp("Poles and Zeros");
disp(p);
disp(z);

%Question 6
IsSysStable = isstable(TransferFunction);
if(IsSysStable == 0)
    disp("System is not BIBO Stable");
else
    disp("System is BIBO Stable");
end

%Question 7
SysWithUnityFeedback = feedback(TransferFunction, 1);


