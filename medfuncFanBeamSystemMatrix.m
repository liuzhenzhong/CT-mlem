function [W_ind,W_dat] = medfuncFanBeamSystemMatrix(theta,N,P_num,delta,SOD,ODD,offset)
% �˺���������������ͶӰ����µ�ϵͳ����
% input parameters:
% theta:ͶӰ�Ƕ�
% N:�����С
% P_num��̽����ͨ������
% delta�������С
% SOD ODD �ֱ�ΪԴ�㵽��ת���ĵľ�������ĵ�̽�����ľ��� offset ƫ��У��ֵ

% output parameters��
% W_ind����ž���  M�У�2*N��
% W_dat�����Ⱦ���  M�У�2*N��

%% -----��Ҫ����ȷ��----- %%
% m0 = 0;
N2 = N^2;
M = length(theta)*P_num; % ͶӰ��������
W_ind = zeros(M,2*N);    % ����
W_dat = zeros(M,2*N);    % ���
t = (-(P_num-1)/2:(P_num-1)/2)*delta;  % ̽���������

% %% -----��ͼ���С��������������ʱ�����Ի���ɨ��ṹͼ����������ͼ----- %%
% if N<=10 && length(theta)<=5
%     x = (-N/2:1:N/2)*delta;
%     y = (-N/2:1:N/2)*delta;
%     plot(x,meshgrid(y,x),'k');
%     hold on;
%     plot(meshgrid(x,y),y,'k')
%     axis([-N/2-5,N/2+5,-N/2-5,N/2+5]);
%     text(0,-0.4*delta,'0');
% end

%% -----ͶӰ����ļ���----- %%
for jj=1:length(theta)
    th = length(theta) - theta(jj) -1;       % ͨ�������������ת���򼴿�
%     th = theta(jj);
    %     sx = SOD*sin(th*pi/180);sy = SOD*cos(th*pi/180);%����Դ����ȷ�� ����ԭ��Ϊ��ת���ģ���ʱ����ת��y�Ὺʼ��ת
    sx = -SOD*sin(th*pi/180);sy = SOD*cos(th*pi/180);
    for ii=1:P_num
        % ���һ������Ȩ���ӵļ���
%         dx = -ODD*sin(th*pi/180) + t(ii)*cos(th*pi/180); %̽Ԫx����ȷ�� ����ԭ��Ϊ��ת���ģ���ʱ����ת��y�Ὺʼ��ת
%         dy = -ODD*cos(th*pi/180) - t(ii)*sin(th*pi /180);%̽Ԫy����ȷ�� ����ԭ��Ϊ��ת���ģ���ʱ����ת��y�Ὺʼ��ת
        dx = ODD*sin(th*pi/180) + (t(ii)+offset)*cos(th*pi/180); %̽Ԫx����ȷ�� ����ԭ��Ϊ��ת���ģ���ʱ����ת��y�Ὺʼ��ת
        dy = -ODD*cos(th*pi/180) + (t(ii)+offset)*sin(th*pi /180);%̽Ԫy����ȷ�� ����ԭ��Ϊ��ת���ģ���ʱ����ת��y�Ὺʼ��ת
        u = zeros(1,2*N);v = zeros(1,2*N);   % ���μ������洢����
        % ����Դ����̽Ԫ����б�ʵ������Ϊ��������������� m��ʾб�� b��ʾ�ؾ�
        % ���Դ����̽Ԫ���ߵ�б��Ϊ�� ƽ����x��
        if abs(sy - dy) < 1e-6   % ƽ����x��
            %             % ������Ӧ������ͼ
            %             if N<=10 && length(theta)<=5
            %                 xx = (-N/2-2:0.01:N/2+2)*delta;
            %                 yy = t(ii);
            %                 plot(xx,yy,'b');
            %                 hold on;
            %             end
            if t(ii)>=N/2*delta||t(ii)<=-N/2*delta         % �����������Χ��ֱ�Ӽ�����һ������
                continue;
            end
            kout = N*ceil(N/2-t(ii)/delta);
            kk = (kout-(N-1)):1:kout;
            u(1:N) = kk;
            v(1:N) = ones(1,N)*delta;
            % ���Դ����̽Ԫ���ߵ�б�ʲ����ڣ���ƽ����Y��
        elseif abs(sx - dx) < 1e-6  % ƽ����Y��
            %             % ������Ӧ������ͼ
            %             if N<=10 && length(theta)<=5
            %                 xx = (-N/2-2:0.01:N/2+2)*delta;
            %                 yy = t(ii);
            %                 plot(xx,yy,'b');
            %                 hold on;
            %             end
            if t(ii)>=N/2*delta||t(ii)<=-N/2*delta         % �����������Χ��ֱ�Ӽ�����һ������
                continue;
            end
            kin = ceil(N/2+t(ii)/delta);
            kk = kin:N:(kin+N*(N-1));
            u(1:N) = kk;
            v(1:N) = ones(1,N)*delta;
            % һ�����,�ȸ��ݵ������ȷ��ֱ�ߵ���������m��b
        else
            %             continue;
            m0 = (sy-dy)/(sx-dx);
            m = abs(m0);
            b = sy - m0*sx;
            %             if m > 0
            y1d = -N/2*delta*m + b;
            y2d = N/2*delta*m + b;
            if (y1d<-N/2*delta && y2d<-N/2*delta)||(y1d>N/2*delta&&y2d>N/2*delta)   % �����������Χ��ֱ�Ӽ�����һ������
                continue;
            end
            %% =====  ȷ�������(xin,yin)�������(xout,yout)������d1 ==== %%
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
                continue;   % ֱ�Ӽ�����һ������
            end
            % ======��������i�������صı�źͳ���====== %
            k = kin;
            c = 0;
            d2 = d1+m*delta;
            while k>=1&&k<=N2              % ����ȥ�������㣬�������±߻�������Χ��ʱ�������¼���
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
                if any(u) == 0 %����������κ�����ֱ�Ӽ�����һ��
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
            %                 if any(u) == 0 %����������κ�����ֱ�Ӽ�����һ��
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
        W_ind((jj-1)*P_num+ii,:) = u;     % ��һ�����ߵ�����������������
        W_dat((jj-1)*P_num+ii,:) = v;     % ��һ�����ߵ�����������������
    end
end
end
