function [ds,dr,dw,da] = TM_dynamics(s,r,w,a,m,I,F,M,Ff,D)
%   检验质量平动和转动动力学

% 输入说明
% s,r,w,a-状态量
% m-TM质量、I-TM惯量
% F-输入控制力、M-输入控制力矩
% Ff-干扰力、D-干扰力矩

% 参数设置
h   = 6671000;                %轨道高度
u_  = 398600e9;               %地心引力常数
w0  = (u_/h^3)^0.5;           %轨道角速度

% 平动系数矩阵
Ar    = [zeros(3) eye(3);
         zeros(3,6)];
Br    = [zeros(3);
         eye(3)/m];
     
% 转动系数矩阵
Ix  = I(1,1);
Iy  = I(2,2);
Iz  = I(3,3);
A41 = -4*w0^2*Ix^-1*(Iy-Iz);
A46 = -w0*Ix^-1*(Iy-Ix-Iz);
A52 = -3*w0^2*Iy^-1*(Ix-Iz);
A63 = -w0^2*Iz^-1*(Iy-Ix);
A64 =  w0*Iz^-1*(Iy-Ix-Iz);
Aw  = [zeros(3) eye(3);
       A41 0 0 0 0 A46;
       0 A52 0 0 0 0;
       0 0 A63 A64 0 0;];
Bw  = [zeros(3);
       1/Ix 0 0;
       0 1/Iy 0;
       0 0 1/Iz];
   
% 平动动力学
x1  = [r;a];
dx1 = Ar*x1+Br*(F+Ff);
dr  = dx1(1:3);
da  = dx1(4:6);

% 转动动力学
x2  = [s;w];
dx2 = Aw*x2+Bw*(M+D);
ds  = dx2(1:3);
dw  = dx2(4:6);

end

