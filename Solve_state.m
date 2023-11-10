function [x,w] = Solve_state(C12,C43,a,b,d0,L)
%   论文中解算位移及转角

%   C12、C43为差分电容
%   d0为初始间距， L为板中心到检验质量中线的距离
%   a为板长，b为板宽
e0 = 8.85e-12;
er = 1;
C0 = er*e0*a*b/d0;
x  = d0*(C12+C43)/(4*C0);
w  = d0*(C43-C12)/(4*C0*L);
end

