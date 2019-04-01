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

C = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1;];

D = [0;0;0;0];

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
disp("Transer Function alpha");
tfAlpha = a(1,1:5);
disp(tfAlpha);
disp(b);
disp("Transer Function alpha dot");
tfAlphaDot = a(2,1:5);
disp(tfAlphaDot);
disp(b);
disp("Transer Function x");
tfx = a(3,1:5);
disp(tfx);
disp(b);
disp("Transer Function x dot");
tfxDot = a(4,1:5);
disp(tfxDot);
disp(b);

%Setting up Transfer Functions
TFSysAlpha = tf(tfAlpha,b);
TFSysAlphaDot = tf(tfAlphaDot,b);
TFSysX = tf(tfx,b);
TFSysXDot = tf(tfxDot,b);
%Question 3

figure;
subplot(2,2,1);
impulse(TFSysAlpha);
subplot(2,2,2);
impulse(TFSysAlphaDot);
subplot(2,2,3);
impulse(TFSysX);
subplot(2,2,4);
impulse(TFSysXDot);

%Question 4
figure;
subplot(2,2,1);
step(TFSysAlpha);
subplot(2,2,2);
step(TFSysAlphaDot);
subplot(2,2,3);
step(TFSysX);
subplot(2,2,4);
step(TFSysXDot);

t = linspace(0,1,1000);
ramp = t;
AlphaRamp = lsim(TFSysAlpha, ramp, t);
AlphaDotRamp = lsim(TFSysAlphaDot, ramp, t);
XRamp = lsim(TFSysX, ramp, t);
XRampDot = lsim(TFSysXDot, ramp, t);
figure;
subplot(2,2,1);
hold on;
plot(t, AlphaRamp, '-');
plot(t, ramp, '--');
xlabel('Time (s)');
ylabel('Output y(t)');
legend('Actual', 'Desired');
title('Alpha Ramp response');

subplot(2,2,2);
hold on;
plot(t, AlphaDotRamp, '-');
plot(t, ramp, '--');
xlabel('Time (s)');
ylabel('Output y(t)');
legend('Actual', 'Desired');
title('Alpha Dot Ramp response');

subplot(2,2,3);
hold on;
plot(t, XRamp, '-');
plot(t, ramp, '--');
xlabel('Time (s)');
ylabel('Output y(t)');
legend('Actual', 'Desired');
title('X Ramp response');

subplot(2,2,4);
hold on;
plot(t, XRampDot, '-');
plot(t, ramp, '--');
xlabel('Time (s)');
ylabel('Output y(t)');
legend('Actual', 'Desired');
title('X Dot Ramp response');

%Question 5 
figure;
subplot(2,2,1);
pzmap(TFSysAlpha);
[pAlpha, zAlpha] = pzmap(TFSysAlpha);
disp("Poles and Zeros Alpha");
disp(pAlpha);
disp(zAlpha);
subplot(2,2,2);
pzmap(TFSysAlphaDot);
[pAlphaDot, zAlphaDot] = pzmap(TFSysAlphaDot);
disp("Poles and Zeros Alpha Dot");
disp(pAlphaDot);
disp(zAlphaDot);
subplot(2,2,3);
pzmap(TFSysX);
[pX, zX] = pzmap(TFSysX);
disp("Poles and Zeros X");
disp(pX);
disp(zX);
subplot(2,2,4);
pzmap(TFSysXDot);
[pXDot, zXDot] = pzmap(TFSysXDot);
disp("Poles and Zeros X Dot");
disp(pXDot);
disp(zXDot);

%Question 6
 IsSysStable = isstable(TFSysAlpha);
 if(IsSysStable == 0)
     disp("Alpha Output is not BIBO Stable");
 else
     disp("Alpha Output is BIBO Stable");
 end
 
 IsSysStable = isstable(TFSysAlphaDot);
 if(IsSysStable == 0)
     disp("Alpha Dot Output is not BIBO Stable");
 else
     disp("Alpha Dot Output is BIBO Stable");
 end
 
 IsSysStable = isstable(TFSysX);
 if(IsSysStable == 0)
     disp("X Output is not BIBO Stable");
 else
     disp("X Output is BIBO Stable");
 end
 
 IsSysStable = isstable(TFSysXDot);
 if(IsSysStable == 0)
     disp("X Dot Output is not BIBO Stable");
 else
     disp("X Dot Output is BIBO Stable");
 end
 
%Question 7

AlphaUnityFeedback = feedback(TFSysAlpha, 1);
AlphaDotUnityFeedback = feedback(TFSysAlphaDot, 1);
XUnityFeedback = feedback(TFSysX, 1);
XDotUnityFeedback = feedback(TFSysXDot, 1);
disp("Alpha Unity");
disp(AlphaUnityFeedback);
disp("Alpha Dot Unity");
disp(AlphaDotUnityFeedback);
disp("X Unity");
disp(XUnityFeedback);
disp("X Dot Unity");
disp(XDotUnityFeedback);
figure;
subplot(2,2,1)
rlocus(AlphaUnityFeedback);
subplot(2,2,2)
rlocus(AlphaDotUnityFeedback);
subplot(2,2,3)
rlocus(XUnityFeedback);
subplot(2,2,4)
rlocus(XDotUnityFeedback);




