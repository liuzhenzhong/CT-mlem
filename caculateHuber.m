function y  = caculateHuber( x1,x2,delta )
%  ����Ϊ���������ֵ�����Ϊ�����������huber����ֵ,deltaΪhuber�����Ľ���ֵ
if abs(x1 - x2) < delta
    y = (x1 - x2);
else
    y = delta * x1;
end
end

