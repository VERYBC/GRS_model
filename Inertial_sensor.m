% 重力参考传感器仿真主程序-生成训练和预测数据
clc
clear
close all

%% 飞行参数
% 传感器参数设置
e0    = 8.85e-12;    % 相对介电常数
er    = 1;           % 真空介电常数
ea    = 1.6e-2;      % 电极板宽度
eb    = 3e-2;        % 电极板长度
L     = 2.15e-2/2;   % 板中心到检验质量中线的距离
tms   = 2.3e-2;      % 检验质量边长的1/2
d     = 4e-3;        % 初始板间距 
m     = 1.96;        % TM质量
I     = [m*tms^2/12 0        0;
         0        m*tms^2/12 0;
         0        0          m*tms^2/12;]; % TM惯量矩阵
     
% 飞行仿真参数设置 
t             = 0;         % 记录时间
dt            = 0.0001;    % 仿真步长
dt_sam        = 1;         % 采样步长
T_control     = 3600;      % 前期稳定时间
T             = 3600*48;   % 稳定后仿真时长
len           = T/dt_sam;  % 采样步数
f1            = 1e-9;      % 干扰力系数
f2            = 1e-5;      % 干扰力矩系数
If_simulation = 0;         % 是否飞行仿真
If_save       = 1;         % 是否保存

%% 故障设置
if_fault      = 1;   % 是否故障:0-无故障、1-偏置故障、2-震荡故障 
fault_C12     = 0;   % 12差分电容故障类型:0-无故障、1-偏置故障、2-震荡故障                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
fault_C34     = 0;   % 34差分电容故障类型
k1            = 0;   % 12差分电容故障系数
k2            = 0;   % 34差分电容故障系数
fault_t_C12   = [];  % 12差分电容故障时间 
fault_t_C34   = [];  % 34差分电容故障时间
fault_len_C12 = [];  % 12差分电容故障时长
fault_len_C34 = [];  % 34差分电容故障时长
fault_dir_C12 = [];  % 12差分电容故障方向
fault_dir_C34 = [];  % 34差分电容故障方向
Add_noise     = 0 ;  % 有色噪声外额外施加噪声

switch if_fault
    case 1 % 偏置故障
        fault_C12     = 1;          % 12差分电容故障类型                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        fault_C34     = 1;          % 34差分电容故障类型
        
        fault_t_C12   = [0.4,0.93];  % 12差分电容故障时间 
        fault_len_C12 = [12,9];      % 12差分电容故障时长
        fault_dir_C12 = [1,-1];      % 12差分电容故障方向
        k1            = [20,10];     % 12差分电容故障系数（5以上）

        
        fault_t_C34   = [0.71];     % 34差分电容故障时间 
        fault_len_C34 = [6];        % 34差分电容故障时长
        fault_dir_C34 = [-1];       % 34差分电容故障方向
        k2            = [15];       % 34差分电容故障系数（5以上）
        
        T             = 1200+300+720; % 预测数据+阈值数据+预测前期数据  
        T_control     = 3600*20;     
        len           = T/dt_sam;
        
        Fault_time = [];
        Fault_this_time    = round([fault_t_C12 fault_t_C34]*len-120); % 记录故障所在时刻，用于替换异常序列
        Fault_next_time    = round([fault_t_C12 fault_t_C34]*len-120) + [fault_len_C12 fault_len_C34]; 
        Fault_this_time    = ceil(Fault_this_time/30);
        Fault_next_time    = ceil(Fault_next_time/30);
        for i = 1:length(Fault_this_time)
            if Fault_next_time(i) ~= Fault_this_time(i)
                Fault_time = [Fault_time Fault_this_time(i) Fault_next_time(i)];
            else
                Fault_time = [Fault_time Fault_this_time(i)];
            end
        end

        Fault_point   = round([fault_t_C12 fault_t_C34]*len-420); % 记录故障具体位置
        Fault_point   = [Fault_point;Fault_point + [fault_len_C12-1 fault_len_C34-1]]; 

    case 2 % 震荡故障
        fault_C12     = 2;        % 12差分电容故障类型                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        fault_C34     = 2;        % 34差分电容故障类型

        fault_t_C12   = [0.33];   % 12差分电容故障时间 （0.26以上）
        fault_len_C12 = [5];      % 12差分电容故障时长
        fault_dir_C12 = [-1];     % 12差分电容故障方向
        k1            = [15];     % 12差分电容故障系数（5以上）

        
        fault_t_C34   = [0.55,0.88];  % 34差分电容故障时间 （0.26以上）
        fault_len_C34 = [10,7];       % 34差分电容故障时长
        fault_dir_C34 = [-1,1];       % 34差分电容故障方向
        k2            = [12,13];      % 34差分电容故障系数（5以上）
        
        T             = 1200+300+720; % 预测数据+阈值数据+预测前期数据  
        T_control     = 3600*21;     
        len           = T/dt_sam;

        Fault_time = [];
        Fault_this_time    = round([fault_t_C12 fault_t_C34]*len-120); % 记录故障所在时刻，用于替换异常序列
        Fault_next_time    = round([fault_t_C12 fault_t_C34]*len-120) + [fault_len_C12 fault_len_C34]; 
        Fault_this_time    = ceil(Fault_this_time/30);
        Fault_next_time    = ceil(Fault_next_time/30);
        for i = 1:length(Fault_this_time)
            if Fault_next_time(i) ~= Fault_this_time(i)
                Fault_time = [Fault_time Fault_this_time(i) Fault_next_time(i)];
            else
                Fault_time = [Fault_time Fault_this_time(i)];
            end
        end

        Fault_point   = round([fault_t_C12 fault_t_C34]*len-420); % 记录故障具体位置
        Fault_point   = [Fault_point;Fault_point + [fault_len_C12 fault_len_C34]]; 

end

%% 加载模型
switch if_fault
    case 0
        load_system('Test_mass.slx')
        set_param('Test_mass','SimulationCommand','start')
        while(string(get_param('Test_mass','SimulationStatus'))== "running")  
              pause(0.5);
        end
    case 1
        load_system('Test_mass_fault.slx')
        set_param('Test_mass_fault','SimulationCommand','start')
        while(string(get_param('Test_mass_fault','SimulationStatus'))== "running")  
              pause(0.5);
        end
    case 2
        load_system('Test_mass_fault.slx')
        set_param('Test_mass_fault','SimulationCommand','start')
        while(string(get_param('Test_mass_fault','SimulationStatus'))== "running")  
              pause(0.5);
        end
end

%% 飞行状态量设置 s-姿态 r-位移 w-角速度 a-加速度
switch If_simulation
    case 1
        s = zeros(len,3);
        r = zeros(len,3);
        w = zeros(len,3);
        a = zeros(len,3);
        F = zeros(len,3); % 控制力
        M = zeros(len,3); % 控制力矩
        % 初值设置
        s(1,:) = [0 0 0];
        r(1,:) = [0 0 0];
        w(1,:) = [0 0 0];
        a(1,:) = [0 0 0];
        % 指令
        sd = [0.1 0.2 0.3];
        rd = [0.1 0.2 0.3];

        %% 飞行仿真（获取姿态及位移） 速度较慢！
        for i=1:len
            Ff = ones(3,1)*cos(0.002*(rand(1)*0.3+0.85)*pi*t)*f1; % 干扰力
            D  = ones(3,1)*cos(0.002*(rand(1)*0.3+0.85)*pi*t)*f2; % 干扰力矩
            % 控制器
            F(i,:) = LQR_control_translation(r(i,:)',a(i,:)',rd');
            M(i,:) = LQR_control_rotation(s(i,:)',w(i,:)',sd',I);
            % 状态更新
            [ds,dr,dw,da] = TM_dynamics(s(i,:)',r(i,:)',w(i,:)',a(i,:)',m,I,F(i,:)',M(i,:)',Ff,D);
            s(i+1,:) = s(i,:) + ds'*dt;
            r(i+1,:) = r(i,:) + dr'*dt;
            w(i+1,:) = w(i,:) + dw'*dt;
            a(i+1,:) = a(i,:) + da'*dt;
            t        = t + dt;
        end
end

%% Simulink读取姿态数据
num = out.s(T_control/dt_sam+2:end,1);
sx  = out.s(T_control/dt_sam+2:end,2);
sy  = out.s(T_control/dt_sam+2:end,3);
sz  = out.s(T_control/dt_sam+2:end,4);
rx  = out.r(T_control/dt_sam+2:end,2);
ry  = out.r(T_control/dt_sam+2:end,3);
rz  = out.r(T_control/dt_sam+2:end,4);

% 传感器状态量空间设置
rx_solve = zeros(len,1);
sz_solve = zeros(len,1);
C        = zeros(len,4);
dC       = zeros(len,2);

% 差分电容传感器噪声
noise1 = Capacitiance_sensor_noise(len,dt_sam,rx(1),1);
noise2 = Capacitiance_sensor_noise(len,dt_sam,rx(1),2);
noise  = [noise1 noise2]; 

% 位移及转角测量值
for i=1:len
    t = t + dt_sam;
    [rx_solve(i),sz_solve(i),C(i,1),C(i,2),C(i,3),C(i,4),dC(i,1),dC(i,2)]...
                                   = Capactance_sensor(rx(i),ry(i),sz(i),...
                                                       ea,eb,d,tms,L,...
                                                       noise(i,:),Add_noise,T,t,...
                                                       fault_C12,fault_t_C12,fault_len_C12,fault_dir_C12,k1,...
                                                       fault_C34,fault_t_C34,fault_len_C34,fault_dir_C34,k2...
                                                       ); 
end

% 标准化
rx_mean = mean(rx_solve);
rx_std  = std(rx_solve);
rx_ = (rx_solve-rx_mean)/rx_std;

%% 保存数据
switch If_save
    case 1
        switch if_fault
            case 0
                rx_save = [num rx_solve];
                col       = {'time','rx_solve'}; 
                result_table=table(rx_save(:,1),rx_save(:,2),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Train_data\rx_solve\rx_solve.csv');

                sz_save = [num sz_solve];
                col       = {'time','sz_solve'}; 
                result_table=table(sz_save(:,1),sz_save(:,2),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Train_data\sz_solve\sz_solve.csv');

                rx_sz_save = [num rx_solve sz_solve];
                col       = {'time','rx_solve','sz_solve'}; 
                result_table=table(rx_sz_save(:,1),rx_sz_save(:,2),rx_sz_save(:,3),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Train_data\rx_sz_solve\rx_sz_solve.csv');

            case 1
                rx_save = [num rx_solve];
                col       = {'time','rx_solve'}; 
                result_table=table(rx_save(:,1),rx_save(:,2),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\偏置故障\rx_solve\rx_solve.csv');

                sz_save = [num sz_solve];
                col       = {'time','sz_solve'}; 
                result_table=table(sz_save(:,1),sz_save(:,2),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\偏置故障\sz_solve\sz_solve.csv');

                rx_sz_save = [num rx_solve sz_solve];
                col       = {'time','rx_solve','sz_solve'}; 
                result_table=table(rx_sz_save(:,1),rx_sz_save(:,2),rx_sz_save(:,3),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\偏置故障\rx_sz_solve\rx_sz_solve.csv');

                fid = fopen('D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\偏置故障\Fault_time.txt','w');
                fprintf(fid,'%g\n',Fault_time);
                fclose(fid);

                fid = fopen('D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\偏置故障\Fault_point.txt','w');
                fprintf(fid,'%g\n',Fault_point);
                fclose(fid);
            case 2
                rx_save = [num rx_solve];
                col       = {'time','rx_solve'}; 
                result_table=table(rx_save(:,1),rx_save(:,2),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\震荡故障\rx_solve\rx_solve.csv');

                sz_save = [num sz_solve];
                col       = {'time','sz_solve'}; 
                result_table=table(sz_save(:,1),sz_save(:,2),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\震荡故障\sz_solve\sz_solve.csv');

                rx_sz_save = [num rx_solve sz_solve];
                col       = {'time','rx_solve','sz_solve'}; 
                result_table=table(rx_sz_save(:,1),rx_sz_save(:,2),rx_sz_save(:,3),'VariableNames',col);
                writetable(result_table, 'D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\震荡故障\rx_sz_solve\rx_sz_solve.csv');

                fid = fopen('D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\震荡故障\Fault_time.txt','w');
                fprintf(fid,'%g\n',Fault_time);
                fclose(fid);

                fid = fopen('D:\Study\Matlab\FDD\Sensor\Inertial sensor\Data\Predict_data\震荡故障\Fault_point.txt','w');
                fprintf(fid,'%g\n',Fault_point);
                fclose(fid);
        end
end


h2=figure(2);
subplot(2,1,1)
f1=plot(rx,'Linewidth',1.5);
hold on
f2=plot(rx_solve,'r--','Linewidth',1.5);
% xticklabel=num2cell(0:T/10:T);
set(gca,'XLim',[0 len]);
% set(gca,'XTicklabel',xticklabel);
set(gca,'child',[f1 f2]);
xlabel('time(s)')
ylabel('d_x(m)')
legend('solve','real')
subplot(2,1,2)
f3=plot(sz,'Linewidth',1.5);
hold on
f4=plot(sz_solve,'r--','Linewidth',1.5);
set(gca,'XLim',[0 len]);
set(gca,'child',[f3 f4])
xlabel('time(s)')
ylabel('\phi_{z}(rad)')
legend('solve','real')
% close(h2)

h3=figure(3);
subplot(2,1,1)
plot(dC(300+121:end,1))
xlabel('time(s)')
ylabel('\DeltaC_1_2(F)')
set(gca,'XLim',[0 len-300-120]);
subplot(2,1,2)
plot(dC(300+121:end,2))
xlabel('time(s)')
ylabel('\DeltaC_3_4(F)')
set(gca,'XLim',[0 len-300-120]);
% close(h3)
