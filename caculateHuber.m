function y  = caculateHuber( x1,x2,delta )
%  输入为两点的像素值，输出为该两点计算后的huber函数值,delta为huber函数的界限值
if abs(x1 - x2) < delta
    y = (x1 - x2);
else
    y = delta * x1;
end
end

