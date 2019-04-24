clear all
clc

load('TDOA_AOA_2LS.mat')
y1 = RMSE;
x1 = x_axis;
figure 
plot(x1,y1,'r-.')
hold on 

load('TDOA_2LS.mat')
y2 = RMSE;
x2 = x_axis;
plot(x2,y2,'b')
legend('TDOA/AOA two-step LS','TDOA two-step LS')
xlabel('probability of location error smaller than ordinate')
ylabel('RMSE Location erro(km)')