%% 240804 truck trailer model
clear
clc
close all
%%
global v tbar L l t0
l = 2.8; % 2.8
L = 5.5; % 5.5
v = - 1; % m/s % -1 ~ -16.2 부터 발산
tbar = 2.0;
t0 = 0.5;
g = 1*t0/pi;


A1 = [-(v*tbar)/(L*t0), 0, 0;
    (v*tbar)/(L*t0), 0, 0;
    (v^2*tbar^2)/(2*L*t0), (v*tbar)/(t0), 0];
A2 = [-(v*tbar)/(L*t0), 0, 0;
    (v*tbar)/(L*t0), 0, 0;
    (g*v^2*tbar^2)/(2*L*t0), (g*v*tbar)/t0, 0];
B1 = [(v*tbar)/(l*t0);
    0;
    0];
B2 = B1;

n_s = size(A1,2);
n_u = size(B1,2);
Q= sdpvar(n_s,n_s,'symmetric');
Y1 = sdpvar(n_u,n_s,'full');
Y2 = sdpvar(n_u,n_s,'full');


LMI1=A1*Q+Q*A1'+B1*Y1+Y1'*B1';

LMI2=A2*Q+Q*A2'+B2*Y2+Y2'*B2';

LMI3=((A1*Q+Q*A1'+B1*Y2+Y2'*B1')/2)+((A2*Q+Q*A2'+B2*Y1+Y1'*B2')/2);

LMI = [LMI1<=0] + [LMI2<=0]+ [LMI3<=0]+ [Q >= 0];

ops=sdpsettings('solver','sedumi','showprogress',0,'verbose',0);
solvesdp(LMI,[],ops)
ch = checkset(LMI) % if eiegenvalues are all negative then ch is all positive

K1 = double(Y1)*inv(double(Q));
K2 = double(Y2)*inv(double(Q));
% save('h_inf_TS_gain_truck.mat','K1','K2')
%%

tf=100;  % final time
ti=0.01;  % runge kutta sample time 
tspan=0:ti:tf;
sample_size = size(tspan,2);
% Assume th \in [-pi pi] 
% th = x(2) + (v*tbar)/(2*L)*x(1);
x(:,1) = [0.5*pi, 0.75*pi,-5]; 
h1 = 1/(1+exp(-2*x(1,1)));
h2 = 1-h1;
U=(h1*K1+h2*K2)*x(:,1);
u_temp(:,1)=U;


cost_w = 0;
    for i=1:sample_size-1
    W(:,i) = 1*exp(-0.1*i*ti)*sin(1*i*ti);
    
    th = x(2,i) + ((v*tbar)/(2*L))*x(1,i);
    if th > pi
        disp('excess pi')
    elseif th < -pi
        disp('excess -pi')
    end
    h1 = 1/(1+exp(-2*x(1,i)));
    h2 = 1-h1;
    
    U=(h1*K1+h2*K2)*x(:,i);
    
    x(:,i+1)=rk5_3(x(:,i),U+W(:,i),ti); 
    u_temp(:,i+1)=U;

    end
figure()
plot(tspan,x(1,:),'r');
hold on
plot(tspan,x(2,:),'b');
plot(tspan,x(3,:),'m');
xlabel('Time')
legend('x_{1}','x_{2}','x_{3}')
grid on

figure()
plot(tspan,u_temp(1,:),'r')
grid on
legend('u_{1}')
xlabel('Time');
% axis([0 10 -10 10])


