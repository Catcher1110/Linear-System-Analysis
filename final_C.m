% This is for Part C
%% This is for (a)
syms w s t tau PhiA(t) A(w);
A(w) = [0,1,0,0;3*w^2, 0,0,2*w; 0,0,0,1; 0,-2*w,0,0];
B = [0 0; 1 0; 0 0; 0 1];
C = [1 0 0 0; 0 0 1 0];
sIA = s*eye(4) - A(w);
% state-transition matrix
PhiA(t) = ilaplace(sIA^-1);
disp(PhiA(t-tau));
% Impulse Response
ImpulseR = C*PhiA(t-tau)*B;
disp(ImpulseR);
% Tranfer Function
TF = C*(sIA^-1)*B;
disp(TF);

%% This is for (b)
syms ts;
Ad = PhiA(ts);
disp(Ad);
Bd = int(PhiA(t),t,0,ts)*B;
disp(Bd);

syms k;
Adk = simplify(Ad^k);
disp(Adk);

%% This is for (d)
Qc = [B, A*B , A^2*B, A^3*B];
disp(rank(Qc));
Qo = [C; C*A; C*A^2; C*A^3];
disp(rank(Qo));
% CEG means controllability gramian matrix
syms theta t0 tf CGM(t0, tf, w)
CGM(t0, tf, w) = simplify(int(PhiA(t0 - tau)*B*transpose(B)*transpose(PhiA(t0 - tau)), tau, t0, tf));
disp(CGM);

%% This is for (e)
syms u(t,t0, tf, theta, w)
x0 = [0;0;theta;0];
xtf = [0;0;(w*tf-theta);0];
u(t,t0, tf, theta, w) = transpose(B)*transpose(PhiA(t0-t))*inv(CGM(t0,tf,w))*(PhiA(t0-tf)*xtf - x0);

syms U1(t, tf) U2(t, tf);
U1(t, tf) = simplify(u(t,0, tf, pi/10, 1/(2*pi)));
U2(t, tf) = simplify(u(t,0, tf, -pi/10, 1/(2*pi)));


%% This is for ploting Energy-tf, it may takes 2-5 minutes for calculation
J1 = zeros(1,191);
J2 = zeros(1,191);
tfx = zeros(1,191);
for tfindex = 0:1:190
    disp(tfindex); % Tell you the progress of the iteration
    tff = 0.1 + 0.01*tfindex;
    tfx(1,tfindex+1) = tff;
    J1(1,tfindex+1) = int(U1(t,tff)' * U1(t, tff), t, 0, tff);
    J2(1,tfindex+1) = int(U2(t,tff)' * U2(t, tff), t, 0, tff);
end

plot(tfx, J1,'r', tfx, J2, 'b');
legend('pi/10','-pi/10')
title('Total Energy Consumption with transfer time')
xlabel('tf /s')
ylabel('J')

%% This is for ploting simulation close-loop system
% tf = 0.1
state_x = zeros(4, 101);
state_x(:,1) = [0;0;pi/10;0];
time_x = zeros(1, 101);
x_dot = zeros(4,1);
tinterval = 0.001;

state_xm = zeros(4, 101);
state_xm(:,1) = [0;0;-pi/10;0];
x_dotm = zeros(4,1);
for time_index = 1:1:100
   disp(time_index);
   ti = tinterval*time_index;
   time_x(1,time_index+1) = ti;
   x_dot = A(1/(2*pi))* state_x(:, time_index) + B*u(ti,0, 0.1, pi/10, 1/(2*pi));
   state_x(:, time_index+1) = state_x(:, time_index) + x_dot * tinterval;
   
   x_dotm = A(1/(2*pi))* state_xm(:, time_index) + B*u(ti,0, 0.1, -pi/10, 1/(2*pi));
   state_xm(:, time_index+1) = state_xm(:, time_index) + x_dotm * tinterval;
end

figure(1);
plot(time_x, state_x(1,:),time_x, state_x(2,:),time_x, state_x(3,:),time_x, state_x(4,:));
title('Simulation of states when tf = 0.1');
legend('x_1','x_2','x_3','x_4');
xlabel('time t /s');

figure(2);
plot(time_x, state_xm(1,:),time_x, state_xm(2,:),time_x, state_xm(3,:),time_x, state_xm(4,:));
title('Simulation of states when tf = 0.1 and \theta = -\pi/10');
legend('x_1','x_2','x_3','x_4');
xlabel('time t /s');

% tf = 1
state_x = zeros(4, 101);
state_x(:,1) = [0;0;pi/10;0];
time_x = zeros(1, 101);
x_dot = zeros(4,1);
tinterval = 0.01;

state_xm = zeros(4, 101);
state_xm(:,1) = [0;0;-pi/10;0];
x_dotm = zeros(4,1);
for time_index = 1:1:100
    disp(time_index);
   ti = tinterval*time_index;
   time_x(1,time_index+1) = ti;
   x_dot = A(1/(2*pi))* state_x(:, time_index) + B*u(ti,0, 1, pi/10, 1/(2*pi));
   state_x(:, time_index+1) = state_x(:, time_index) + x_dot * tinterval;
   
   x_dotm = A(1/(2*pi))* state_xm(:, time_index) + B*u(ti,0, 1, -pi/10, 1/(2*pi));
   state_xm(:, time_index+1) = state_xm(:, time_index) + x_dotm * tinterval;
end

figure(1);
plot(time_x, state_x(1,:),time_x, state_x(2,:),time_x, state_x(3,:),time_x, state_x(4,:));
title('Simulation of states when tf = 1');
legend('x_1','x_2','x_3','x_4');
xlabel('time t /s');

figure(2);
plot(time_x, state_xm(1,:),time_x, state_xm(2,:),time_x, state_xm(3,:),time_x, state_xm(4,:));
title('Simulation of states when tf = 1 and \theta = -\pi/10');
legend('x_1','x_2','x_3','x_4');
xlabel('time t /s');
% tf = 2
state_x = zeros(4, 201);
state_x(:,1) = [0;0;pi/10;0];
time_x = zeros(1, 201);
x_dot = zeros(4,1);
tinterval = 0.01;

state_xm = zeros(4, 201);
state_xm(:,1) = [0;0;-pi/10;0];
x_dotm = zeros(4,1);
for time_index = 1:1:200
    disp(time_index);
   ti = tinterval*time_index;
   time_x(1,time_index+1) = ti;
   x_dot = A(1/(2*pi))* state_x(:, time_index) + B*u(ti,0, 2, pi/10, 1/(2*pi));
   state_x(:, time_index+1) = state_x(:, time_index) + x_dot * tinterval;
   
   x_dotm = A(1/(2*pi))* state_xm(:, time_index) + B*u(ti,0, 2, -pi/10, 1/(2*pi));
   state_xm(:, time_index+1) = state_xm(:, time_index) + x_dotm * tinterval;
end

figure(1);
plot(time_x, state_x(1,:),time_x, state_x(2,:),time_x, state_x(3,:),time_x, state_x(4,:));
title('Simulation of states when tf = 2');
legend('x_1','x_2','x_3','x_4');
xlabel('time t /s');

figure(2);
plot(time_x, state_xm(1,:),time_x, state_xm(2,:),time_x, state_xm(3,:),time_x, state_xm(4,:));
title('Simulation of states when tf = 2 and \theta = -\pi/10');
legend('x_1','x_2','x_3','x_4');
xlabel('time t /s');

%% This is for (f)
F = [-0.5 0 0 0; 
    0 -2.5 0 0;
    0 0 -0.5+1i 0;
    0 0 0 -0.5-1i];
Lhat = [1 1 1 1;
        1 1 1 1];
disp(rank([Lhat; Lhat*F; Lhat*F^2; Lhat*F^3]));

T = lyap3(transpose(A), -F, -transpose(C) * Lhat);
L = transpose(Lhat * inv(T));
disp(L);
disp(eig(simplify(A-L*C)));

%% This is for g
F = [-1 0 0 0; 
    0 -2 0 0;
    0 0 -3 0;
    0 0 0 -4];
Lhat = [1 1 1 1;
        1 1 1 1];
disp(rank([Lhat; Lhat*F; Lhat*F^2; Lhat*F^3]));
T = lyap3(transpose(A), -F, -transpose(C) * Lhat);
H1 = transpose(Lhat * inv(T));
disp(simplify(H1));
disp(eig(simplify(A-H1*C)));

F = [-7 0 0 0; 
    0 -8 0 0;
    0 0 -9 0;
    0 0 0 -10];
Lhat = [1 1 1 1;
        1 1 1 1];
disp(rank([Lhat; Lhat*F; Lhat*F^2; Lhat*F^3]));
T = lyap3(transpose(A), -F, -transpose(C) * Lhat);
H2 = transpose(Lhat * inv(T));
disp(simplify(H2));
disp(eig(simplify(A-H2*C)));

T = [eye(4);eye(4)];
F1 = A - H1*C;
F2 = A - H2*C;
Fhat = [F1, zeros(4,4);
    zeros(4,4),F2];
disp(vpa(simplify(det([T, expm(Fhat*2)*T]))));



