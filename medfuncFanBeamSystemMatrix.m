function [W_ind,W_dat] = medfuncFanBeamSystemMatrix(theta,N,P_num,delta,SOD,ODD,offset)
% 此函数用来计算扇束投影情况下的系统矩阵
% input parameters:
% theta:投影角度
% N:矩阵大小
% P_num；探测器通道个数
% delta；网格大小
% SOD ODD 分别为源点到旋转中心的距离和中心到探测器的距离 offset 偏移校正值

% output parameters：
% W_ind：编号矩阵  M行，2*N列
% W_dat：长度矩阵  M行，2*N列

%% -----重要参数确定----- %%
% m0 = 0;
N2 = N^2;
M = length(theta)*P_num; % 投影的总条数
W_ind = zeros(M,2*N);    % 长度
W_dat = zeros(M,2*N);    % 序号
t = (-(P_num-1)/2:(P_num-1)/2)*delta;  % 探测器坐标号

% %% -----当图像较小，射线条数较少时，可以画出扫描结构图，这是网格图----- %%
% if N<=10 && length(theta)<=5
%     x = (-N/2:1:N/2)*delta;
%     y = (-N/2:1:N/2)*delta;
%     plot(x,meshgrid(y,x),'k');
%     hold on;
%     plot(meshgrid(x,y),y,'k')
%     axis([-N/2-5,N/2+5,-N/2-5,N/2+5]);
%     text(0,-0.4*delta,'0');
% end

%% -----投影矩阵的计算----- %%
for jj=1:length(theta)
    th = length(theta) - theta(jj) -1;       % 通过这个来控制旋转方向即可
%     th = theta(jj);
    %     sx = SOD*sin(th*pi/180);sy = SOD*cos(th*pi/180);%射线源坐标确定 坐标原点为旋转中心，逆时针旋转，y轴开始旋转
    sx = -SOD*sin(th*pi/180);sy = SOD*cos(th*pi/180);
    for ii=1:P_num
        % 完成一条射线权因子的计算
%         dx = -ODD*sin(th*pi/180) + t(ii)*cos(th*pi/180); %探元x坐标确定 坐标原点为旋转中心，逆时针旋转，y轴开始旋转
%         dy = -ODD*cos(th*pi/180) - t(ii)*sin(th*pi /180);%探元y坐标确定 坐标原点为旋转中心，逆时针旋转，y轴开始旋转
        dx = ODD*sin(th*pi/180) + (t(ii)+offset)*cos(th*pi/180); %探元x坐标确定 坐标原点为旋转中心，逆时针旋转，y轴开始旋转
        dy = -ODD*cos(th*pi/180) + (t(ii)+offset)*sin(th*pi /180);%探元y坐标确定 坐标原点为旋转中心，逆时针旋转，y轴开始旋转
        u = zeros(1,2*N);v = zeros(1,2*N);   % 单次计算结果存储矩阵
        % 根据源点与探元连线斜率的情况分为三种情况进行讨论 m表示斜率 b表示截距
        % 如果源点与探元连线的斜率为零 平行于x轴
        if abs(sy - dy) < 1e-6   % 平行于x轴
            %             % 画出对应的射线图
            %             if N<=10 && length(theta)<=5
            %                 xx = (-N/2-2:0.01:N/2+2)*delta;
            %                 yy = t(ii);
            %                 plot(xx,yy,'b');
            %                 hold on;
            %             end
            if t(ii)>=N/2*delta||t(ii)<=-N/2*delta         % 如果超出网格范围，直接计算下一条射线
                continue;
            end
            kout = N*ceil(N/2-t(ii)/delta);
            kk = (kout-(N-1)):1:kout;
            u(1:N) = kk;
            v(1:N) = ones(1,N)*delta;
            % 如果源点与探元连线的斜率不存在，即平行于Y轴
        elseif abs(sx - dx) < 1e-6  % 平行于Y轴
            %             % 画出对应的射线图
            %             if N<=10 && length(theta)<=5
            %                 xx = (-N/2-2:0.01:N/2+2)*delta;
            %                 yy = t(ii);
            %                 plot(xx,yy,'b');
            %                 hold on;
            %             end
            if t(ii)>=N/2*delta||t(ii)<=-N/2*delta         % 如果超出网格范围，直接计算下一条射线
                continue;
            end
            kin = ceil(N/2+t(ii)/delta);
            kk = kin:N:(kin+N*(N-1));
            u(1:N) = kk;
            v(1:N) = ones(1,N)*delta;
            % 一般情况,先根据点的坐标确定直线的两个参数m和b
        else
            %             continue;
            m0 = (sy-dy)/(sx-dx);
            m = abs(m0);
            b = sy - m0*sx;
            %             if m > 0
            y1d = -N/2*delta*m + b;
            y2d = N/2*delta*m + b;
            if (y1d<-N/2*delta && y2d<-N/2*delta)||(y1d>N/2*delta&&y2d>N/2*delta)   % 如果超出网格范围，直接计算下一条射线
                continue;
            end
            %% =====  确认入射点(xin,yin)、出射点(xout,yout)及参数d1 ==== %%
            if y1d<=N/2*delta&&y1d>=-N/2*delta&&y2d>N/2*delta
                yin = y1d;
                d1 = yin-floor(yin/delta)*delta;
                kin = N*floor(N/2-yin/delta)+1;
                yout = N/2*delta;
                xout = (yout-b)/m;
                kout = ceil(xout/delta)+N/2;
                
            elseif y1d<=N/2*delta && y1d>=-N/2*delta && y2d>=-N/2*delta && y2d<N/2*delta
                yin = y1d;
                d1 = yin-floor(yin/delta)*delta;
                kin = N*floor(N/2-yin/delta)+1;
                yout = y2d;
                kout = N*floor(N/2-yout/delta)+N;
                
            elseif y1d<-N/2*delta&&y2d>N/2*delta
                yin = -N/2*delta;
                xin = (yin-b)/m;
                d1 = N/2*delta+(floor(xin/delta)*delta*m+b);
                kin = N*(N-1)+N/2+ceil(xin/delta);
                yout = N/2*delta;
                xout = (yout-b)/m;
                kout = ceil(xout/delta)+N/2;
                
            elseif y1d<-N/2*delta && y2d>=-N/2*delta && y2d<N/2*delta
                yin = -N/2*delta;
                xin = (yin-b)/m;
                d1 = N/2*delta+(floor(xin/delta)*delta*m+b);
                kin = N*(N-1)+N/2+ceil(xin/delta);
                yout = y2d;
                kout = N*floor(N/2-yout/delta)+N;
            else
                continue;   % 直接计算下一条射线
            end
            % ======计算射线i穿过像素的编号和长度====== %
            k = kin;
            c = 0;
            d2 = d1+m*delta;
            while k>=1&&k<=N2              % 像素去遍历计算，当像素下边还在网格范围内时进行如下计算
                c = c+1;
                if d1 >= 0&&d2 > delta
                    u(c)=k;
                    v(c)=(delta-d1)*sqrt(m^2+1)/m;
                    if k>N&&k~=kout
                        k=k-N;
                        d1=d1-delta;
                        d2 = d1+m*delta;
                    else
                        break;
                    end
                elseif d1>=0&&d2==delta
                    u(c) = k;
                    v(c) = delta*sqrt(m^2+1);
                    if k>N&&k~=kout
                        k = k-N+1;
                        d1 = 0;
                        d2 = d1+m*delta;
                    else
                        break;
                    end
                elseif d1>=0&&d2<delta
                    u(c) = k;
                    v(c) = delta*sqrt(m^2+1);
                    if k~=kout
                        k=k+1;
                        d1=d2;
                        d2 = d1+m*delta;
                    else
                        break;
                    end
                elseif d1<=0&&d2>=0&&d2<=delta
                    u(c) = k;
                    v(c) = d2*sqrt(m^2+1)/m;
                    if k~=kout
                        k=k+1;
                        d1=d2;
                        d2=d1+m*delta;
                    else
                        break;
                    end
                elseif d1<=0&&d2>delta
                    u(c)=k;
                    v(c) = delta*sqrt(m^2+1)/m;
                    if k>N&&k~=kout
                        k=k-N;
                        d1=-delta+d1;
                        d2 = d1+m*delta;
                    else
                        break;
                    end
                end
            end
            if m0 < 0
                u_temp = zeros(1,2*N);
                if any(u) == 0 %如果不经过任何网格，直接计算下一条
                    continue;
                end
                ind = u>0;
                for k=1:length(u(ind))
                    r = rem(u(k),N);
                    if r==0
                        u_temp(k) = u(k)-N+1;
                    else
                        u_temp(k) = u(k)-2*r+N+1;
                    end
                end
                u = u_temp;
            end
            %         elseif m<0
            %                 u_temp = zeros(1,2*N);
            %                 if any(u) == 0 %如果不经过任何网格，直接计算下一条
            %                     continue;
            %                 end
            %                 ind = u>0;
            %                 for k=1:length(u(ind))
            %                     r = rem(u(k),N);
            %                     if r==0
            %                         u_temp(k) = u(k)-N+1;
            %                     else
            %                         u_temp(k) = u(k)-2*r+N+1;
            %                     end
            %                 end
            %                 u = u_temp;
            %                   continue;
            %             end
        end
        W_ind((jj-1)*P_num+ii,:) = u;     % 将一条射线的情况存入输出矩阵中
        W_dat((jj-1)*P_num+ii,:) = v;     % 将一条射线的情况存入输出矩阵中
    end
end
end
