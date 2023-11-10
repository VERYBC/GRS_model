function C = Nonparallel_capacitance(d0,a,w,b)
% 非平行平板电容
% d0为两平板中心距，a为板长，b为板宽，w为转角
e0 = 8.85e-12;
er = 1;
C  = b*e0*er*log(1+(a*tan(w)/(d0-a*tan(w)/2)))/tan(w);
end

