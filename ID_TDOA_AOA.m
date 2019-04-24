%%
%本程序为TDOA/AOA综合定位结果,,使用解法为Two-step LS
% BS为一条直线，斜率为1，共4个BS
% TDOA噪声标准差为0.1km，角度为1度
% MS位置（0.25，0.75）
%%
clear all
clc
%% 初始化
% BS
x1 = [0,0];
x2 = [5,5];
x3 = [10,10];
x4 = [-5,-5];
BS = [x1;x2;x3;x4];

% MS
x_real = 0.25;
y_real = 0.75;

c = 3*10e5;
sigma_t = 0.1/c;       %可调参数，c*sigma_t为TDOA标准差
s_beta = 1*(pi/180);

M = 3;
%噪声协方差矩阵
Q = (0.5*eye(M)+0.5*ones(M))*(sigma_t^2);

N = 10000;
postion = zeros(2,N);
for num = 1:N
    n_TDOA = (c*sigma_t).*randn(M,1);
    n_beta = s_beta*randn(1,1);
    %% 搭建测量模型
    % 定义1基站为home station，MS到1距离记作R1
    R1 = sqrt((x_real - x1(1))^2+(y_real - x1(2))^2);

    % 产生TDOA测量值，第i个基站到MS距离与R1的差
    R = zeros(M,1);
    for i =1:M
        R(i) = sqrt((x_real - BS(i+1,1))^2+(y_real - BS(i+1,2))^2)-R1+n_TDOA(i);
    end
    % 产生AOA测量值
    beta = atan((y_real - x1(2))/(x_real - x1(1)))+n_beta;
    
    G = zeros(M,2);
    for i = 1:M
        G(i,1) = BS(i+1,1);
        G(i,2) = R(i);
    end
    G = (-1).*G;

    h = zeros(M,1);
    K = zeros(M,1);
    for i = 1:M
        K(i+1) = BS(i+1,1)^2 + BS(i+1,2)^2;
        h(i) = R(i)^2-K(i+1);
        h(i) = 0.5*h(i);
    end
    %% 根据测量值估值
    GQ = G'*inv(Q);
    z = inv(GQ*G)*GQ*h;
    
    w = z(1);
    r1 =z(2);
    py = (w+sqrt(2*r1^2-w^2))/2;
    px = w-py;
    %估计B，再迭代一次

    B = zeros(M,M);
    for i =1:M
        temp = sqrt((px - BS(i+1,1))^2+(py - BS(i+1,2))^2)-R1;
        B(i,i) = temp;
    end
%     D1 = sqrt((px - BS(1,1))^2+(py - BS(1,2))^2);

    fai = c^2.*B*Q*B;
    GF = G'*inv(fai);
    z_2 = inv(GF*G)*GF*h;
    w = z_2(1);
    r1 =z_2(2);
    y = (w+sqrt(2*r1^2-w^2))/2;
    x = w-y;
    postion(1,num)=x;
    postion(2,num)=y;
end
% N次定位结果的定位误差距离
erro = zeros(N,1);
for i = 1:N
    erro(i) = sqrt((postion(1,i)-x_real)^2+(postion(2,i)-y_real)^2);
end

%% 作图
% 改变误差阈值，统计阈值内的点的RMSE
erro_max = max(erro);
count = 0;   %统计误差阈值内定位结果个数
p = 0:0.1:ceil(erro_max); %距离误差阈值 ，单位km
for i = 1:length(p)
    if i ==1
        x_axis(i)=0;
        RMSE(i)=0;
    else  
        temp = 0;
        count = 0;
        for j =1:N
            if erro(j)<=p(i)
                count = count+1; 
                temp = temp+erro(j)^2;
            end
        end
        x_axis(i) = count/N;
        RMSE(i) =  sqrt(temp/count);
    end
end

save('1-D_TDOA.mat','x_axis','RMSE')
figure 
plot(x_axis,RMSE)
title('1-D Array BS TDOA')
xlabel('probability of location error smaller than ordinate')
ylabel('RMSE Location erro(km)')
hold on
%% 加入AOA信息
TA_postion = zeros(2,N);
Gl = zeros(3,2);
Gl(1:2,1:2)=eye(2);


for i = 1:N
    n_beta = s_beta*randn(1,1);
    % 产生AOA测量值
    beta = atan((y_real - x1(2))/(x_real - x1(1)))+n_beta;
    Gl(3,1) = -sin(beta);
    Gl(3,2) = cos(beta);
    ml = [postion(1,i);postion(2,i);0];
    D1 = sqrt((postion(1,i) - BS(1,1))^2+(postion(2,i) - BS(1,2))^2);
    zp = inv(Gl'*Gl)*Gl'*ml;
    TA_postion(1,i)=zp(1);
    TA_postion(2,i)=zp(2);
end

% N次定位结果的定位误差距离
TA_erro = zeros(N,1);
for i = 1:N
    TA_erro(i) = sqrt((TA_postion(1,i)-x_real)^2+(TA_postion(2,i)-y_real)^2);
end

% 改变误差阈值，统计阈值内的点的RMSE
TA_erro_max = max(TA_erro);
count = 0;   %统计误差阈值内定位结果个数
p = 0:0.1:ceil(TA_erro_max); %距离误差阈值 ，单位km
for i = 1:length(p)
    if i ==1
        TA_x_axis(i)=0;
        RMSE(i)=0;
    else  
        temp = 0;
        count = 0;
        for j =1:N
            if TA_erro(j)<=p(i)
                count = count+1; 
                temp = temp+TA_erro(j)^2;
            end
        end
        TA_x_axis(i) = count/N;
        TA_RMSE(i) =  sqrt(temp/count);
    end
end


plot(TA_x_axis,TA_RMSE,'-.')
title('1-D Array BS TDOA/AOA')
xlabel('probability of location error smaller than ordinate')
ylabel('RMSE Location erro(km)')
legend('TDOA two-step LS','TDOA/AOA two-step LS')

