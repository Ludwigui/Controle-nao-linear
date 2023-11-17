clc
clearvars
yalmip('clear');
format long e
%% Sistema de controle de nível em um reservatório cônico.
A2 = 7.9*10^(-3); %m^2
g = 9.8; %m/s^2

hbar = 5;
%Sistema original:
% dot(h) = b*u/(h^2)+a*sqrt(h)/(h^2)
b = 4/(pi*0.6^2);
a = 4*A2*sqrt(2*g)/(pi*0.6^2);

%% Definindo valor de controle no equilíbrio:

dbar = a/b*sqrt(hbar);
%% É necessário transladar o sistema para o ponto.

%x1 = h - hbar e u(t)=d(t)-dbar
% Sendo hdot = b/(h^2)*d(t)-a*sqrt(h)/(h^2) ----- d(t) a abertura da
% válvula.
%d(t)=u(t)+dbar

%% Portanto:
%dotx1 = -a/(sqrt((x1+hbar)^3)) + b*dbar/((x1+hbar)^2) + b*u(t)/((x1+hbar)^2)
%% Adiciona-se um integrador no sistema:
% dotx2 = x1
%% Definição do sistema
%Abordagem LPV:
%(a0 - theta1*a1 + theta2*a2)*x1 = -a*x1/(sqrt((x1+hbar)^3)*x1) + b*dbar*x1/((x1+hbar)^2*x1)
% Ou seja, a0 = 0, a1 = -a, a2 = b*dbar

%(b0 + theta3*b1)*u(t) = b*u(t)/((x1+hbar)^2)
% Ou seja, b0 = 0 e b1 = b

%Definição dos thetas:
%theta1 = 1/(sqrt((x1+hbar)^3)*x1)        -200.4 <= theta1 <= 0.0063
%theta2 = 1/((x1+hbar)^2*x1)              -2004 <= theta2 <= 0.002
%theta3 = 1/((x1+hbar)^2)                 0.01 <= theta3 <= 10000  

PS = [-200.4 -2004 0.01;
      -200.4 0.002 0.01;
      -200.4 -2004 10000;
      -200.4 0.002 10000;
       0.0063 -2004 0.01;
       0.0063 0.002 0.01;
       0.0063 -2004 10000;
       0.0063 0.002 10000];

a0 = 0;
a1 = -a;
a2 = b*dbar;

b0 = 0;
b1 = b;

Amatriz = [0 0;
           1 0];
Bmatriz = [0;0];


Ad = Amatriz;
Bd = Bmatriz;
mod = 4.99;
limX1 = (1/mod)*[1;0];
limX2 = (1/mod)*[0;1];
%% System dimensions
nx = size(Amatriz,1);
nu = size(Bmatriz,2);
[nl,nv] = size(PS);
%% Guaranteed decay rate
r = 1;
%% Decision variables
Q = sdpvar(nx,nx);
W0 = sdpvar(nu,nx,'full');
W1 = sdpvar(nu,nx);
W2 = sdpvar(nu,nx);
mu = sdpvar(1,1);
options = sdpsettings('verbose',0,'solver','sdpt3');
%% LMIs definition
 LMIs = [Q >= 0] + [1 - limX1'*Q*limX1 >= 0] + [1 - limX2'*Q*limX2 >= 0];
 for i=1:nv
     Ad(1,1) = Ad(1,1) - a1*PS(i,1) + a2*PS(i,2);
     W = W0 + PS(i,1)*W1 + PS(i,2)*W2;
     Bd(1,1) = Bd(1,1) + b1*PS(i,3);
     lmi = Ad*Q + Bd*W + r*Q;
     LMIs = LMIs + [(lmi+lmi') <= 0];
 end
%% --------------------------------------------------------------------------
% Calling the solver
result = optimize(LMIs,-geomean(Q),options);
%--------------------------------------------------------------------------
% Recovering the results
Q = double(Q);
P = inv(Q);
W0 = double(W0);
W1 = double(W1);
W2 = double(W2);
K0 = W0*inv(Q);
K1 = W1*inv(Q);
K2 = W2*inv(Q);
%--------------------------------------------------------------------------
% Result validation
test_LMI = check(LMIs);
disp(result);
disp(test_LMI);
disp('The stabilization with ROA estimate is');
display(Q);
display(K0);
display(K1);
display(K2);
save('s8_ex2_solution','P','Q','K0','K1','K2','r');


