clc
clearvars
clear all
close all
yalmip('clear');
%% Definicao dos parâmetros e limites do sistema

a = 0.86;
b = 3.53;
hbar = 6;
d = 0.6;

x1a = linspace(-4,4,100);
fx1 = (-a)./(sqrt((x1a+hbar).^3)) + b*d./((x1a+hbar).^2);
fx2 = b./((x1a+hbar).^2);

figure
plot(x1a,fx1)
grid on
title('Análise dos limites de f_{xa}')

figure
plot(x1a,fx2)
grid on
title('Análise dos limites de f_{xb}')

%% Definindo o sistema
% fx1_a = (a0+a1*theta1)x1
a0 = -0.0288;
a1 = -0.0262;
% fx1_b = (b0+b1*theta1)*u(t)
b0 = 0.4680;
b1 = -0.4320;
% Inicialmente, defino a matriz aumentada sem a0, a1, b0 e b1
A = [0 0; 1 0];
B = [0;0];
V = [-1 -1 1 1;
    1 -1 1 -1]; %Vértices do politopo
limX1 = 4; % -4 <= X1 <= 4
alfa1 = [1/limX1; 0];

%% Taxa de decaimento
r = 0.025;
%% Dimensões 
nx = size(A,1);
nu = size(B,2);
nv = size(V,2);
%% Variáveis de decisão
Q = sdpvar(nx,nx);
W = sdpvar(nu,nx);
M = [1 2.5*W;
    2.5*W' Q];
options = sdpsettings('verbose',0,'solver','sdpt3');
%% Definição das LMIs
LMIs = (Q>=0) + (1-alfa1'*Q*alfa1>=0);
for i=1:nv
    A(1,1) = a0+a1*V(1,i);
    B(1,1) = b0+b1*V(2,i);
    lmi = A*Q+B*W+r*Q;
    LMIs = LMIs + ((lmi+lmi')<=0) + (M>=0);
end
%% Calling the solver
result = optimize(LMIs,-geomean(Q),options); 
%% Recovering the results
Q = double(Q);
P = inv(Q);
W = double(W);
K = W*inv(Q);
M = double(M);
%% Result validation
test_LMI = check(LMIs);
disp(result);
disp(test_LMI);
disp('The stabilization with ROA estimate is');
display(Q);
display(K);
%% Pegando os ganhos para simulação
k1 = K(1,1);
k2 = K(1,2);
%% Condições iniciais
x1ini=0.05;
x2ini=-0.5;
TStop = 150;
sim('simuNoLinearControl')
% Plotando a região de atração do sistema
% ROA estimation
figure
title('ROA')
bq=16;
[x1s,x2s] = meshgrid(-bq:0.01:bq,-bq:0.01:bq);
[nx1,mx1]=size(x1s);
[nx2,mx2]=size(x2s);
for i=1:nx1
    for j=1:nx2
        X = [x1s(i,j) x2s(i,j)]';
        z(i,j) = X'*P*X;
    end
end
hold on
contour(x1s,x2s,z,[1 1],'b','LineWidth',2)
hold on
plot(data(:,2),data(:,3),'r','LineWidth',2)
hold off
grid on
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
% Ploting the real domain of attraction
% It is needed to make the time simulation in advance
%hold on
%plot(x.signals.values(:,1),x.signals.values(:,2))

%% Condições iniciais

figure
plot(data(:,1),data(:,2))
hold on
plot(data(:,1),data(:,3))
legend('x_1','x_2')
xlabel('Tempo [s]')
grid on

figure
plot(datau(:,1),datau(:,2))
legend('u(t)=d(t)-\bar(d)')
xlabel('Tempo [s]')
grid on
ylim([-0.25 0.05])







