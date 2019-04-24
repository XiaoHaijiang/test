%%
clear all
clc
%% 
c = 3e8;
% anchor���ʼ��
x1 = [0,0];
x2 = [5,0];
x3 = [5,5];
x4 = [0,5];
BS = [x1;x2;x3;x4];
BS = BS.*1e3;
M =4;

%agent���ʼ����agent��ֻ��1������
x_real = 3e3;
y_real = 3e3;          %�ο���λ��

% ���������ź�
D = zeros(M,1);      %M��anchor�㵽agent��ľ���
tao = zeros(M,1);     %M��anchor�㵽agent�㴫��ʱ��
for i =1:M
      D(i) = sqrt((x_real - BS(i,1))^2+(y_real - BS(i,2))^2);
      tao(i) = D(i)/c;   
end

fc = 400e6;       %��Ƶ400MHz
fs = 2*fc;           %����Ƶ��
f = 20e3;          %ԭʼ�ź�Ƶ��20k
t = 0:1/fs:1/(2*f)-1/fs;

r_signal = zeros(M,length(t));

A = 1;
SNR = 10 ;            %����ȵ�λdB
sigma = sqrt(0.5*A^2/(10^(SNR/10)));

NUM = 1000;            %��������ֵ����
noise = sigma*randn(NUM,length(t));

R = zeros(M-1,NUM);         %TDOA��ʵֵ
T = zeros(M-1,NUM);         %TDOA����ֵ
s_beta = 1*(pi/180);
n_beta = s_beta*randn(NUM,1);
beta = zeros(NUM,1);        %AOA����ֵ

% ��ͨ�˲���
wp=0.1*pi;    
ws=0.3*pi;
deltaw=ws-wp;
N=ceil(6.6*pi/deltaw);
wdhamm=hamming(N)';
wc=(wp+ws)/2;      %����Ƶ��Ϊfs/10 = fc/5 = 80MHz
hd=ideallp(wc,N);
h=hd.*wdhamm;
fai = [0,0,-pi/2,0];
for num =  1:NUM; 
    % ����AOA����ֵ
    beta(num) = atan((y_real/1e3 - x1(2))/(x_real/1e3 - x1(1)))+n_beta(num);
    for i =1:M
      ti = t - tao(i).*ones(1,length(t));
      s0 = abs(cos(2*pi*f.*ti));
      r_signal(i,:) = A*(s0.*cos(2*pi*fc.*ti+fai(i)))+noise(num,:);
    %   r_signal(i,:) = s0;
    end
    %% ���
    H = fft(h,length(t));
    s = zeros(M,length(t));
    for i =1:M
       rs = r_signal(i,:).*cos(2*pi*fc*t);
       Rs = fft(rs,length(t));
       S = Rs.*H;
       s(i,:) = abs(real(ifft(S)));
    end
    % ���ʱ����pi��λģ����ԭ����cos(2*pi*fc*tao),��˲��ð벨��Ϊ�����ź�
    
    %% ��ش����õ�TDOA����ֵ
    for i =1:M-1
       R(i,num) = tao(i+1)-tao(1);
       T(i,num) = tdoa_measure(s(i+1,:),s(1,:),fs);
    end
    
end

%% ��λ�㷨
model.BS = BS/1e3;
model.sigma = sigma;
model.s_beta = s_beta;
model.f = f;
p_tdoa = zeros(2,NUM);
p_tdoa_aoa = zeros(2,NUM);
for i = 1:NUM
  p_tdoa(:,i) = TDOA(model,T(:,i));
  p_tdoa_aoa(:,i) = TDOA_AOA(model,T(:,i),beta(i));
end

real(1)=x_real/1e3;
real(2)=y_real/1e3;
[x_tdoa,RMSE_toda] = cal_erro(p_tdoa,real);
[x_tdoa_aoa,RMSE_toda_aoa] = cal_erro(p_tdoa_aoa,real);
figure 
plot(x_tdoa_aoa,RMSE_toda_aoa,'r-.')
hold on 
plot(x_tdoa,RMSE_toda,'b')

legend('TDOA/AOA two-step LS','TDOA two-step LS')
xlabel('probability of location error smaller than ordinate')
ylabel('RMSE Location erro(km)')




