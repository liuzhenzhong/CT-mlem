function F = medfuncMedianMap( W_ind,W_dat,N,F0,P,irt_num ,beta)
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   ���������
%   ���������
N2 = N^2;
M = length(P); % ������������
F = F0;
k=0;
temp = zeros(N2,1);
while(k<irt_num)
    disp(['��' num2str(irt_num) '�ε�������' num2str(k+1) '�ε���']);
    %% ==��ʵ��ͶӰֵ�͹���ͶӰֵ�ı�ֵ����=== %%
    image = reshape(F,N,N);
    U = getHuber(image,0.05);
    U = reshape(U,N2,1);
    %     U = caculateUx(image);
    %     U = reshape(U,N2,1);
    %     image = medfilt2(image,[3,3]);
    %     image = reshape(image,N2,1);
    R = zeros(M,1);
    delta = zeros(M,1);
    for ii=1:M
        % ������߲�ͨ���κ����أ���������
        if any(W_ind(ii,:)) == 0
            continue;
        end
        w = zeros(1,N2);%ϵͳ�����һ��������
        for jj=1:2*N
            m=W_ind(ii,jj);
            if m>0&&m<=N2
                w(m)=W_dat(ii,jj);
            end
        end
        proj = w*F;%ǰ��ͶӰ
        delta(ii) = P(ii) - proj;  % ͶӰ��ǰ��ͶӰ�Ĳ�ֵ����
        R(ii)=P(ii)./proj;%��ֵ������һ������һ����ֵ
    end
    % ====== ���ֵ������ȫ1�����ķ�ͶӰ ========= %
    back_proj = zeros(N2,1);back_proj0=zeros(N2,1);
    for ii=1:M %���еı��������з�ͶӰ����
        label = W_ind(ii,:);
        data = W_dat(ii,:);
        if any(label)~=0
            ind = label>0;
            index = label(ind);
            back_proj(index) = back_proj(index)+R(ii)*data(ind).';%��ֵ�����ķ�ͶӰ
            back_proj0(index) = back_proj0(index)+data(ind).';%+ beta*U(index);%ȫ1�����ķ�ͶӰ
        end
    end
    % =========ͼ�����һ��============== %
    ind = back_proj0>0;
    %     F(ind) = F(ind).*(back_proj(ind)./back_proj0(ind));   % ��ͨ��ml_em�㷨
    %     F(ind) = F(ind).*(back_proj(ind)./(back_proj0(ind)+beta*U(ind)));    % ��������Ʒ����������֪ʶ���ؽ�
%     F(ind) = F(ind).*(back_proj(ind)./(back_proj0(ind)+beta*((F(ind)-img(ind))./img(ind))+beta*U(ind)));% ������ֵ������ؽ�
    F(ind) = F(ind).*(back_proj(ind)./(back_proj0(ind)+50*temp(ind)+beta*U(ind)));
    k=k+1;   
    % ����betaֵ  ����Ӧ����betaֵ
    %     son = sqrt(dot(delta,delta));
    %     mother = sqrt(dot(U,U));
    %     beta = log((son/(mother+0.1))+1);
    % ����ÿ�ε����Ľ��
    result_save = reshape(F,N,N)';
    writebin(['median_' num2str(k) '.bin'],result_save);  
    %% ���һ�ε�������һ���˲����� 
%     img0 = MedianFilter(result_save,N,100,4000);
%     img = reshape(img0,N2,1);
    img = medfilt1(F,16);
    temp = (F-img)./(img+0.001);
end
end

