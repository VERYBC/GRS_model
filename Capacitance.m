function C= capacitance(x,y,w,d0,s,L,a,b)
%   计算电容（注意敏感轴正向）
%   来源：ICCES会议论文

%   x为敏感轴方向位移，y为非敏感轴方向位移，w为敏感轴方向转角
%   d0为初始间距，s为检验质量边长的一半， L为板中心到检验质量中线的距离
%   a为板长，b为板宽
%   中心距d1
d1 = d0+s-x-s*sec(w)-y*tan(w)+L*tan(w);
C  = Nonparallel_capacitance(d1,a,w,b);
end
