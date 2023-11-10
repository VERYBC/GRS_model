 function [x_solve,w_solve,C1,C2,C3,C4,C12,C43] = Capactance_sensor(x,y,w,...
                                                                   a,b,d0,s,L,...
                                                                   noise,Add_noise,T,t,...
                                                                   fault_C12,fault_t_C12,fault_len_C12,fault_dir_C12,k1,...
                                                                   fault_C34,fault_t_C34,fault_len_C34,fault_dir_C34,k2...
                                                                   )
%   电容传感器测量值_设置了故障
%   输入状态，输出测量值
%   0为无故障；1为偏置故障；2为震荡故障；3；
%   k为故障系数大小
%----C1
C1 = Capacitance(x,y,w,d0,s,L,a,b);
%----C2
C2 = Capacitance(-x,y,-w,d0,s,L,a,b);
%----C3
C3 = Capacitance(-x,y,w,d0,s,L,a,b);
%----C4
C4 = Capacitance(x,y,-w,d0,s,L,a,b);

% 计算差分电容
C12 = C1 - C2 + noise(1) + normrnd(0,Add_noise*5e-18);
C43 = C4 - C3 + noise(2) + normrnd(0,Add_noise*5e-18);

% 故障设置
switch fault_C12
    case 1
        for i = 1:length(fault_t_C12)
            if t>=round(T*fault_t_C12(i)) && t<round(T*fault_t_C12(i)+fault_len_C12(i))
                C12 = C12  + 2e-17*k1(i)*fault_dir_C12(i);
            end
        end
    case 2
        for i = 1:length(fault_t_C12)
            if t>=round(T*fault_t_C12(i)) && t<round(T*fault_t_C12(i)+fault_len_C12(i))
                C12 = C12  + normrnd(0,2e-17) * k1(i);
            end
        end
end
switch fault_C34
    case 1
        for i = 1:length(fault_t_C34)
            if t>=round(T*fault_t_C34(i)) && t<round(T*fault_t_C34(i)+fault_len_C34(i))
                C43 = C43  + 2e-17*k2(i)*fault_dir_C34(i);
            end
        end
    case 2
        for i = 1:length(fault_t_C34)
            if t>=round(T*fault_t_C34(i)) && t<round(T*fault_t_C34(i)+fault_len_C34(i))
                C43 = C43  + normrnd(0,2e-17)* k2(i);
            end
        end    
end

% 解算状态量
[x_solve,w_solve] = Solve_state(C12,C43,a,b,d0,L);

end

