% EES-40 2023 - Controle Moderno Lab4 
% Controlador utilizando placa USB NI-6009 e MATLAB R2016a
% Script para projetar controle discreto do kit: DLQR com filtro de Kalman (KF) = DLQG  
% controle estocástico DLQG com realimentação de estimativa do estado AUMENTADO 
% estado aumentado com  x_i: ação integral para rastrear referência degrau em malha fechada 
% pseudomedida zk=x_i usada no KF
% realização no espaço de estado.
% rodar após o script de projeto ter rodado 
clear all; 
clc;
load Kdlqr % DLQR gain
load Ts    % sampling time [s]
load Gaugssd % zoh discretization of augmented plant model
Aaugd=Gaugssd.A;
Baugd=Gaugssd.B(:,1);
Gaugd=Gaugssd.B(:,2);
Caugd=Gaugssd.C(1,:);
C_xid=Gaugssd.C(2,:);

disp('Início');
daqreset;
DeviceID = 'Dev2'; % Confirmar número do Dev no software da National

%% Setup:
% Tempos:
T=Ts; % intervalo de amostragem em segundos adotado no projeto do controlador
Tempo_Experiencia = 30; % Tempo total da experiência em segundos

% Ganho do AmpOp
Ganho_Fisico = 5; % Hardware comercial configurado na placa borne do AmpOp 
                  % os valores têm tolerância na fabricação

% Inicialização:
t = 0; % Tempo inicial
time = [];
u = [];     % DLQR control
y = [];     % actual plant output
r = [];     % reference step from signal generator
z = [];     % pseudomeasurement z=xik=integral(r-y)=integral(ek)
inov = [];   % innovation z-zhat
Sinov = []; % KF-computed inovation variance 

P0=diag([.01 1 1 .01]);  % initial augmented state estimation error covariance matrix 
x0=zeros(4,1);           % initial augmented state estimation error mean
uk=0;                    % initial control sample
zk=0;                    % initial pseudomeasurement sample

Pprop=zeros(4,4);       % memory space for propagated augmented covariance matrix
Pupdt=P0;               % updated augmented covariance matrix
Gain=zeros(4,1);        % memory space for Kalman gain vector
xhat_prop=zeros(4,1);   % memory space for filter propagated augmented state vector
xhat_updt=zeros(4,1);   % memory space for filter updated augmented state vector

qd=.25*Ts;           % discrete zero-mean Gaussian model noise wd(k) variance is PSD*Ts 
qdf=qd;              % filter tuning parameter
rd=.0001/Ts;         % discrete zero-mean Gaussian pseudomeasurement noise v(k) variance is PSD/Ts
rdf=rd;              % filter tuning parameter

% Placa NI-6009:
s = daq.createSession('ni');
s.Rate = 1/T;
addAnalogInputChannel(s,DeviceID,0:1,'Voltage'); % 2 entradas analógicas do A/D
addAnalogOutputChannel(s,DeviceID,0:1,'Voltage'); % 2 saídas analógicas do D/A
%% Loop de controle:
tic; % start
while t < Tempo_Experiencia 
    clc
    disp(sprintf('Controlador funcionando em tempo real. Aguarde por %u segundos.',Tempo_Experiencia - t));
    
    % KF covariance propagation, innovation variance, and gain
    Pprop=Aaugd*Pupdt*Aaugd'+Baugd*qdf*Baugd';  
    Sinovk=C_xid*Pprop*C_xid'+rdf;              
    Gain=Pprop*C_xid'/Sinovk;  
    
    % Amostrando os sinais na entrada do conversor A/D
    In = inputSingleScan(s); 
    yk = In(:,1); % Sinal de saída da planta = ai0
    ref = In(:,2); % Sinal de referência para o sistema = ai1
    
    ek=ref-yk;     % output tracking error
    zk=zk+Ts*ek;   % updated zoh discretization of integral of ek = pseudomeasurement
    
    % KF mean propagation
    xhat_prop=Aaugd*xhat_updt+Baugd*uk+Gaugd*ref;
    % KF pseudomeasurement innovation - should be white sequence
    inovk=zk-C_xid*xhat_prop;   
    % KF mean update
    xhat_updt=xhat_prop+Gain*inovk;
    
    % DLQG stochastic control
    uk=-Kdlqr*xhat_updt;
    
    % Níveis de saturação do sinal de controle em -15V e 15V
    uk = min(uk,15);
    uk = max(-15,uk);
  
    % Adequação do sinal de saída do DAC limitado de 0 a 5 V que deve ser
    % utilizado com o AmpOp na configuração diferencial com ganho externo
    % Vo = (R2/R1)*(Vb - Va), onde Va = a0 e Vb = a1
    if uk > 0 % Sinal de saída do AmpOp deve ser positivo
        a1 = uk/Ganho_Fisico; a0 = 0;
       else % Sinal de saída do AmpOp deve ser negativo
        a1 = 0; a0 = -uk/Ganho_Fisico;
    end  % entrada para planta fica no intervalo -3V a 3V
    
    % Escrevendo o sinal de controle nos conversores D/A
    outputSingleScan(s,[a0,a1]);
    
    % KF covariance update
    Pupdt=(eye(4)-Gain*C_xid)*Pprop;  
   
   % Atualizando o tempo
    t = t + T;
  
    % Armazenando valores 
    u = [u;uk];
    y = [y;yk]; % plant output
    r = [r;ref];
    time = [time;t];
    z = [z;zk]; % pseudomeasurement
    inov=[inov;inovk]; % pseudomeasurement innovation
    Sinov=[Sinov;Sinovk]; % KF-computed pseudomeasurement variance
   
    % Sincronizando o intervalo de amostragem
    pause(T-toc);
    tic;
end
% Zerando as saídas do conversor D/A:
outputSingleScan(s,[0,0]);
disp('Acabou');

% Liberando o AD/DA
release(s);

%% Desenhando gráficos
save inov.mat inov
save Sinov.mat Sinov
figure(1);
plot(time,y,'-r',time,r,'-.k');
title('Sinais de referência e saída do servo');
grid; xlabel('t(s)'); ylabel('y(V) e r(V)');

figure(2);
stairs(time,u,'-b')
title('Sinal de controle do servo');
grid; xlabel('kT(s)'); ylabel('u(V)');

% Innovation: just one realization and KF-computed stdev
aa4=Sinov;
dd4=inov;
figure(3);
plot(time,sqrt(Sinov),'x',time,-sqrt(Sinov),'x');
hold;plot(time,inov,'o');    
title('innovation sequence (should be white) - KF stdev(+/-) and a typical realization');
figure(4);
[inovautocorr,lags]=xcorr(inov,'normalized'); %Matlab v2016a: usar 'coeff'
stem(lags,inovautocorr);
title('innovation sequence (should be white) - a typical realization - autocorrelation');


% END



