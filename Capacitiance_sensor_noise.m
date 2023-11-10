function [noise] = Capacitiance_sensor_noise(length,dt_sam,delta,i)
%   生成电容传感器噪声（时间序列）
%   来源【Capacitive sensing of test mass motion with nanometer precision over millimeter-wide sensing gaps for space-borne
%         gravitational reference sensors】

% 输入
% dt_sam采样步长,delta位移,第i个电容传感器(x:1-2; y:3-4; z:5-6)

% 功率谱密度PSD
A   = 0.78e-18; % 低频噪声
B   = (31+rand(1)*4-2)*1e-6*delta*1e-6;
switch i 
    case 1
        C = (0.57+(rand(1)*2-1)*0.03)*1e-18;
    case 2
        C = (0.54+(rand(1)*2-1)*0.05)*1e-18;
    case 3
        C = (0.85+(rand(1)*2-1)*0.03)*1e-18;
    case 4
        C = (0.49+(rand(1)*2-1)*0.08)*1e-18;
    case 5
        C = (0.31+(rand(1)*2-1)*0.03)*1e-18;
    case 6
        C = (0.64+(rand(1)*2-1)*0.04)*1e-18;        
end

fs  = 1000;     % 噪声频率
f   = 1/dt_sam; % 采样频率
fk  = linspace(0.0001,10,length*(fs/f)/2);
PSD = A^2+B^2*(1e-3*(fk.^-1))+C^2*(1e-3*(fk.^-1));
% 生成单边频谱幅值
PSD_double = PSD/2;
% 生成频谱幅值
Ax = sqrt(PSD_double*(length*(fs/f)/2)*fs);
% 生成频谱
fik = pi*randn(1,length*(fs/f)/2);   
xk  = Ax.*exp(1i*fik); 
xk  = [xk flip(xk)];
% 生成数据
xm = ifft(xk);                      
noise_1000hz = real(xm);
noise = zeros(length,1);
for i=1:length
    noise(i) = noise_1000hz(i*(fs/f));
end

