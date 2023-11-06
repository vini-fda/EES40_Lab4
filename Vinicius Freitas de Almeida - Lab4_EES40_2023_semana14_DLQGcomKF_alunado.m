% EES-40 2023-2 Lab4 semana 14 Jacques Waldmann novembro 2023
% Controle ótimo DLQG de sistema SISO 
% Polos de malha fechada do sistema contínuo selecionados via root locus simétrico (LGR-S) 
% Robustez do controle ótimo contínuo LQR - função de sensibilidade e margens de estabilidade
% Verificação da banda passante da malha fechada para discretização e projeto discreto DLQR 
% Estimador: filtro de Kalman => sistema variante no tempo 
% Desempenho DLQG e o compromisso com o esforço de controle e a passagem de ruído de medida 

% Realização em cascata de primeira ordem com  modal de segunda ordem
% G(s)= 10/(s+10) * 100/(s+0,5+10j)/(s+0,5-10j) = Y(s)/U(s) é a malha aberta
clear all;clc;close all
% cada grupo com realização {A,B,C} distinta
A=[-10 0 0;-10 -.5 10;0 -10 -.5];
B=[10;0;0];
C=[0 0 1];
Gss=ss(A,B,C,0);
[Gtfnum Gtfden]=ss2tf(A,B,C,0);
Gtf=tf(Gtfnum,Gtfden); % malha aberta, sem ação integral
step(Gtf)
%% augmented state xaug dynamics 
% includes xi: time integral of r-y (output error wrt to reference) 
% xaug_dot=Aaug*xaug+Baug*u+Gaug*r
% y=Caug*xaug - plant model output measurement
Aaug=[A zeros(3,1);-C 0];
Baug=[B;0];
Gaug=[zeros(3,1); 1]; % reference r input matrix  
Caug=[C 0]; % physical output y
% C_xi is the infinite-horizon time integral of the output tracking error 
% for use in the quadratic cost function J_LQR
C_xi=[0 0 0 1]; % pseudomeasurement y'=x_i
if rank(ctrb(Aaug,Baug))==length(Aaug) 
    X=['sistema aumentado contínuo é controlavel'];
    disp(X);
else disp('sistema continuo nao é controlavel');
end
s=tf('s'); % Laplace transform - variable s for later use in transfer functions  
%% lugar geométrico das raízes simétrico (LGR-S)
% symmetric root locus (SRL) LGR-S to select closed-loop poles according to rho^(-1)
% minimize J_LQR=integral (xi^2 + rho u^2)dt from t=0 to +inf
[numAug,denAug]=ss2tf(Aaug,Baug,C_xi,0);
Gaugtf=tf(numAug,denAug); % malha aberta, inclui ação integral, com "saída" xi ponderada em J_LQR
% lqr function yields Klqr(1:4)=[K  -Ki] such that (be aware of the negative sign)
%      u= - Klqr * [x'  xi']'= - Klqr * xaug is the control signal optimizing J_LQR 
% while driving the augmented linear dynamics and tracking the reference.  
% LQR design assumes the actual state is available 
% Closed-loop pole selection according to symmetric root locus with output xi in J_LQR  
num1=numAug; % numerador de Gaugtf(-s)
den1=denAug.*[1 -1 1 -1 1]; % denominador de Gaugtf(-s), que tem xi como "saída"
num_s=conv(numAug,num1); % numerador de Gaugtf(s)*Gaugtf(-s)
den_s=conv(denAug,den1); % denominador de Gaugtf(s)*Gaugtf(-s)
sys_s=tf(num_s,den_s); % open loop tf Gaugtf(s)*Gaugtf(-s) for SRL
rlocus(sys_s); % choose with mouse acceptable gain from the stable closed-loop poles in SRL asymptotes
GainLGR_S=input(' selected SRL gain (inverse of rho seen in J-cost functional ');
rho=1/GainLGR_S
Klqr=lqr(Aaug,Baug,C_xi'*C_xi,rho); % lqr com peso em xi e rho=-1/GainLGR_S no controle, sendo GainLGR_S obtido com LGR simétrico para atender requisitos
%% closing the loop with continuous-time LQR optimal state feedback gain
Ki=-Klqr(4);
K=Klqr(1:3);
Aaugmf=[Aaug-Baug*Klqr]; % closed-loop dynamics matrix
X=['state feedback gain vector K ',num2str(K(1:3)),' output tracking error integral gain Ki ',num2str(Ki)];
disp(X);
% closed-loop measurement and pseudomeasurement output signals in LQR feedback design
% read bandwidth wb from closed-loop (Y(jw)/R(jw) plot at the top 
SysMF=ss(Aaugmf,Gaug,[Caug;C_xi],zeros(2,1)); 
% closed-loop, augmented-state LQR design
% notice: input is reference r, output signals are y and xi=xaug(4)
damp(SysMF);
p=eig(Aaugmf);
% read bandwidth wb [rd/s] in closed-loop Y(jw)/R(jw) Bode magnitude plot
figure;bodemag(SysMF); 
%% model-based steady state analysis: closed-loop response in steady state 
step_mag=5; % reference step magnitude [V] 
X=['ref magnitude ',num2str(step_mag),' [V]'];disp(X);
MM=[A B;-C 0];
xureg=inv(MM)*[zeros(3,1);-step_mag] % plant model steady state is xureg(1:3) and control xureg(4)
xireg=1/Ki*(xureg(4)+K*xureg(1:3)) % augmented state component xi in steady state
%% sampling time choice and ZOH plant discretization for DLQR design
wb=input('-3 dB bandwidth frequency [rd/s] from Bode closed loop Y(j\omega)/R(j\omega) ');
fb=wb/2/pi; % Hz
% NI-board maximum D/A sampling frequency is 150Hz
Ts=1/(20*fb); 
X=['suggested sampling time Ts ',num2str(Ts),' [s]'];
disp(X);
Ts=input(' enter suitable sampling time [s]');
save Ts.mat Ts % for use in NI board driver script

disp('desired z-plane closed loop poles ')
Closed_z=exp(p*Ts)
% augmented plant model, state space realization
Gaugss=ss(Aaug,[Baug Gaug],[Caug;C_xi],zeros(2,2));
Gaugssd=c2d(Gaugss,Ts,'zoh');  % ZOH discretization
save Gaugssd.mat Gaugssd
% zoh discretization of augmented plant model
Aaugd=Gaugssd.A;
Baugd=Gaugssd.B(:,1);
Gaugd=Gaugssd.B(:,2);
Caugd=Gaugssd.C(1,:);
C_xid=Gaugssd.C(2,:);
if rank(ctrb(Aaugd,Baugd))==length(Aaugd) 
    X=['sistema discreto controlavel'];
    disp(X);
else disp('sistema discreto nao controlavel');
end
%% DLQR design
% discrete-time state feedback design optimizes cost J_DLQR=Sum(x_i[k]^2 + rho*u[k]^2) with u[k] = -Kdlqr*xaug[k]
[Kdlqr,~,Eig_d]=dlqr(Aaugd,Baugd,C_xid'*C_xid,rho); % dlqr pondera x_i e rho=-1/GainLGR_S no controle, sendo GainLGR_S obtido com LGR simétrico para atender requisitos
save Kdlqr.mat Kdlqr;
disp('dlqr closed-loop poles');
Eig_d

X=['discrete-time state feedback gain vector Kdlqr ',num2str(Kdlqr)];
disp(X);

Amfd=Aaugd-Baugd*Kdlqr;
disp('Amfd=Aaugd-Baugd*Kdlqr closed-loop poles - compare with Eig_d');
pd=eig(Amfd) % polos de malha fechada
%% discrete-time closed-loop simulation of DLQR design
% reference is r=step_mag*1(t) and output {y;xi]
SysMFd=ss(Amfd,Gaugd,[Caug;C_xi],zeros(2,1),Ts);
step_mag=5; 
opt = stepDataOptions('StepAmplitude',step_mag);
tfinal=5; % final time tf=5[s]

% xstep is the timeline of xaug; ystep is the timeline of xi=xaug(4)
[ystep,tstep,xstep]=step(SysMFd,tfinal,opt);  
sizek=length(tstep); % length to be used in Kalman filter steps

figure; % sinal xaug1(t) em malha fechada
stem(tstep,xstep(:,1));
title(['xaug(1) em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

figure; % sinal xaug2(t) em malha fechada
stem(tstep,xstep(:,2));
title(['xaug(2) em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

figure; % sinal xaug(3) em malha fechada
stem(tstep,xstep(:,3));
title(['xaug(3) em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

figure; % sinal xaug(4) em malha fechada
stem(tstep,xstep(:,4));
title(['xaug(4) em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

figure; % output signals: y(k);xi(k)  
stem(tstep,ystep(:,1:2));
title(['saída y e pseudomedida x_I em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V.s]']);
legend('y','x_I');
xlabel('segundos');
grid;

figure;   % sinal de controle em malha fechada 
xstep_transp=xstep'; % linha temporal do estado aumentado em malha fechada é aqui transposto
cntrl=-Kdlqr*xstep_transp;
stem(tstep,cntrl);
title(['controle DLQR [V] em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

%% DLQR design stability margins with augmented state xaug feedback is dictated by 
% the open-loop transfer function tfL(z), z=exp(jw.Ts)
% u=Kdlqr*xaug as output and negative feedback is assumed
% projeto DLQR,assegura estabilidade exponencial assegurada.
% inspecionar quão próximo do ponto crítico passa o modelo nominal 
SysMAd=ss(Aaugd,Baugd,Kdlqr,0,Ts);
figure;
nyquistplot(SysMAd); % use zoom in critical point in the Nyquist plot and see margins
title('Nyquist plot - open loop tfL(z)[dB] at z=exp(j*w*Ts)');
figure;
margin(SysMAd);
legend('SysMAd(z), z=exp(jw.Ts): open-loop DLQR'); 
% para verificar se projeto é ou não é robusto
% ver se tfL(z),z=exp(jw.Ts) se aproxima muito do ponto crítico - pode se tornar instável com modelo incerto e envolver ponto crítico 
[NumD,DenD]=ss2tf(Aaugd,Baugd,Kdlqr,0); 
tfLd=tf(NumD,DenD,Ts); % open loop discrete tf with u(k)=Kdlqr*xaug(k) as output and negative feedback is assumed 
Dlqr=1+tfLd; % distance from critical point in complex plane (-1+0j) to open loop tfL(z),z=exp(jw.Ts)
Slqr=1/Dlqr; % sensitivity tf
figure;
bodemag(Slqr);
title('Magnitude da função de sensibilidade |S_{dlqr}(exp(Ts.j\omega))|[dB]'); 
figure;
bodemag(Dlqr);
title('função distância ao ponto crítico |D_{dlqr}(exp(Ts.j\omega))|[dB]'); 
%% observability test 
if rank(obsv(Aaugd,C_xi))==length(Aaugd) 
    X=['sistema discreto observavel'];
    disp(X);
else disp('sistema discreto nao observavel');
end
%% Lab 4 - Malha fechada com estimação via filtro de Kalman
% KF here estimates the entire augmented state vector
% Closing the loop with the Kalman filter - covariances and innovation variance analysis 
%                                            and tuning the filter noise parameters qdf, rdf
x0=[0 0 0 0]';       % ground-truth plant initial state                        
P0=diag([.01 1 1 .01]);  % initial state estimation error covariance matrix 

qd=.25*Ts;           % discrete zero-mean Gaussian model noise wd(k) variance is PSD*Ts 
qdf=qd;              % filter tuning parameter
rd=.00005/Ts;        % .0001/Ts discrete zero-mean Gaussian pseudomeasurement noise v(k) variance is PSD/Ts
rdf=rd;              % filter tuning parameter
N=sizek;             % N here is filter steps length - integer
M=1;%M=1000;         % M here is number of independent realizations in the Monte Carlo simulation

% memory space for logging various arrays
Pprop=zeros(4,4,N);       % memory space for propagated covariance matrix
Pupdt=zeros(4,4,N);       % memory space for updated covariance matrix
Gain=zeros(4,N);        % memory space for Kalman gain vector
x=zeros(4,N);           % memory space for ground-truth state vector
z=zeros(1,N);           % memory space for pseudomeasurement z=y´=x_i vector
xhat_prop=zeros(4,N);   % memory space for propagated filter state vector
xhat_updt=zeros(4,N);   % memory space for updated filter state vector

inov=zeros(1,N,M);      % memory space for inovation for Monte Carlo statistics
uhat=zeros(1,N,M);      % memory space for control for MC statistics 
oute=zeros(1,N,M);      % memory space for output tracking error e=stepmag-Caugd*x

Sinov=zeros(1,N);       % memory space for KF inovation variance 

% KF covariances and innovation variance
Pprop(:,:,1)=Aaugd*P0*Aaugd'+Baugd*qdf*Baugd';    % P(1|0) uses P0,qdf
Sinov(1,1)=C_xid*Pprop(:,:,1)*C_xid'+rdf;         % S(1|0) uses rdf     
Gain(:,1)=Pprop(:,:,1)*C_xid'/Sinov(1,1);         % L(1)
Pupdt(:,:,1)=(eye(4)-Gain(:,1)*C_xid)*Pprop(:,:,1);    % P(1|1)
 
 % operação iterativa - KF covariances and innovation variance
 for k=2:N
  Pprop(:,:,k)=Aaugd*Pupdt(:,:,k-1)*Aaugd'+Baugd*qdf*Baugd';     % P(k|k-1)
  Sinov(1,k)=C_xid*Pprop(:,:,k-1)*C_xid'+rdf;                % S(k|k-1)      
  Gain(:,k)=Pprop(:,:,k)*C_xid'/Sinov(1,k);             % L(k)
  Pupdt(:,:,k)=(eye(4)-Gain(:,k)*C_xid)*Pprop(:,:,k);   % P(k|k)
 end  % end of KF covariance and innovation variance loop

 % KF-computed variances -  plots
steps=1:N;

% KF xaug(1) stdev - propagation and update
% KF xaug(1) estimation gain
a1(1:N)=Pprop(1,1,1:N);
b1(1:N)=Pupdt(1,1,1:N);
figure;
plot(tstep,sqrt(a1),'x',tstep,-sqrt(a1),'x');hold;
plot(tstep,sqrt(b1),'o',tstep,-sqrt(b1),'o');
title('xaug(1) estimation uncertainty - prop and updt KF-computed filter stdev');
figure;
c1(1:N)=Gain(1,1:N);
stem(tstep,c1,'x');title('KF-computed xaug(1) estimation Gain');

% KF xaug(2) stdev - propagation and update
% KF xaug(2) estimation gain
a2(1:N)=Pprop(2,2,1:N);
b2(1:N)=Pupdt(2,2,1:N);
figure;
plot(tstep,sqrt(a2),'x',tstep,-sqrt(a2),'x');hold;    
plot(tstep,sqrt(b2),'o',tstep,-sqrt(b2),'o');
title('xaug(2) estimation uncertainty - prop and updt KF-computed stdev');
figure;
c2(1:N)=Gain(2,1:N);
stem(tstep,c2,'x');title('KF-computed xaug(2) estimation gain');

% KF-computed xaug(3) uncertainty stdev - propagation and update
% KF xaug(3) estimation gain
a3(1:N)=Pprop(3,3,1:N);
b3(1:N)=Pupdt(3,3,1:N);
figure;
plot(tstep,sqrt(a3),'x',tstep,-sqrt(a3),'x');hold;    
plot(tstep,sqrt(b3),'o',tstep,-sqrt(b3),'o');  
title('xaug(3) estimation uncertainty - prop and updt KF-computed stdev');
figure;
c3(1:N)=Gain(3,1:N);
stem(tstep,c3,'x');title('KF xaug(3) estimation gain');

% KF-computed xaug(4) uncertainty stdev - propagation and update
% KF xaug(4) estimation gain
a4(1:N)=Pprop(4,4,1:N);
b4(1:N)=Pupdt(4,4,1:N);
figure;
plot(tstep,sqrt(a4),'x',tstep,-sqrt(a4),'x');hold;    
plot(tstep,sqrt(b4),'o',tstep,-sqrt(b4),'o');  
title('xaug(4) estimation uncertainty - prop and updt KF-computed stdev');
figure;
c4(1:N)=Gain(4,1:N);
stem(tstep,c4,'x');title('KF xaug(4) estimation gain');
%% KF 
% reference step from signal generator is r=step_mag*1(t)
for j=1:M        % Monte Carlo independent j-th realization loop
 % first iteration:
 % ground-truth process and pseudomeasurement, tracking error and control
 oute0=step_mag-Caugd*x0;            % initial plant output tracking error 
 x0hat=mvnrnd(x0,P0)';    % initial filter state estimate 4x1 drawn from multivariable Gaussian density N(x0,P0)
 u0hat=-Kdlqr*x0hat;                  % initial control signal uses initial filter estimate 
 % true stochastic state uses initial true state(=0)
 x(:,1)=Aaugd*x0+Baugd*u0hat+Gaugd*step_mag+Baugd*sqrt(qd)*randn;  
 z(1,1)=C_xid*x(:,1)+sqrt(rd)*randn;                % true stochastic x_i pseudomeasurement
 oute(1,1,j)=step_mag-Caugd*x(:,1);                  % tracking error 
 % KF
 xhat_prop(:,1)=Aaugd*x0hat+Baugd*uhat(1,1,j)+Gaugd*step_mag;
 xhat_updt(:,1)=xhat_prop(:,1)+Gain(:,1)*(z(1,1)-C_xid*xhat_prop(:,1));
 inov(1,1,j)=z(1,1)-C_xid*xhat_prop(:,1);          % j-th realization of initial inovation
 uhat(1,1,j)=-Kdlqr*xhat_updt(:,1);             % control signal uses filter estimate

 % KF operação iterativa - ground-truth e filtro
 for k=2:N
  % ground-truth process and pseudomeasurement
  x(:,k)=Aaugd*x(:,k-1)+Baugd*uhat(1,k-1,j)+Gaugd*step_mag+Baugd*sqrt(qd)*randn;   % true state propagates
  oute(1,k,j)=step_mag-Caugd*x(:,k);  % tracking error by the end of true state propagation with control signal
  z(1,k)=C_xid*x(:,k-1)+sqrt(rd)*randn;  % true pseudomeasurement 
  % KF 
  xhat_prop(:,k)=Aaugd*xhat_updt(:,k-1)+Baugd*uhat(1,k-1,j)+Gaugd*step_mag;
  xhat_updt(:,k)=xhat_prop(:,k)+Gain(:,k)*(z(1,k)-C_xid*xhat_prop(:,k));
  xerr_prop(:,k,j)=xhat_prop(:,k)-x(:,k);     % j-th realization of k-th-step propagated state estimation error
  xerr_updt(:,k,j)=xhat_updt(:,k)-x(:,k);     % j-th realization of k-th-step updated state estimation error 
  inov(1,k,j)=z(1,k)-C_xid*xhat_prop(:,k);    % j-th realization of k-th-step inovation
  uhat(1,k,j)=-Kdlqr*xhat_updt(:,k);  % control signal uses filter estimate and reference
 end  % end of KF step loop
end  % end of Monte Carlo independent realization loop
%% gráficos 
figure;stem(tstep,-oute+step_mag);title('Plant model output - one realization');
figure;stem(tstep,uhat);title('Control - one realization');
%% Innovation: just one realization and KF-computed stdev
aa4(1:N)=Sinov(1,1:N);
dd4(1:N)=inov(1,:,1);
figure;
plot(tstep,sqrt(aa4),'x',tstep,-sqrt(aa4),'x');
hold;plot(tstep,dd4,'o');    
title('innovation - KF stdev(+/-) and one realization');
figure;
[inovautocorr,lags]=xcorr(dd4,'normalized');
stem(lags,inovautocorr);
title('innovation - one realization - autocorrelation');

