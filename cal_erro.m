function [x_axis,RMSE] = cal_erro(postion,real)
    N = size(postion,2);
    x_real = real(1);
    y_real = real(2);
    % N�ζ�λ����Ķ�λ������
    erro = zeros(N,1);
    for i = 1:N
        erro(i) = sqrt((postion(1,i)-x_real)^2+(postion(2,i)-y_real)^2);
    end
    % ��ͼ
    % �ı������ֵ��ͳ����ֵ�ڵĵ��RMSE
    erro_max = max(erro);
    count = 0;   %ͳ�������ֵ�ڶ�λ�������
    p = 0:0.1:floor(erro_max); %���������ֵ ����λkm
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
            if count==0
                x_axis(i) = 0;
                RMSE(i) = 0;
            else
                x_axis(i) = count/N;
                RMSE(i) =  sqrt(temp/count);
            end
        end
    end