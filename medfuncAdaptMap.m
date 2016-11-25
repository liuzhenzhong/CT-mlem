function F = medfuncAdaptMap( W_ind,W_dat,N,F0,P,irt_num ,beta)
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   ���������
%   ���������
N2 = N^2;
M = length(P); % ������������
F = F0;
k=0;
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
    F(ind) = F(ind).*(back_proj(ind)./(back_proj0(ind)+beta*U(ind)));    % ��������Ʒ����������֪ʶ���ؽ�
%     F(ind) = F(ind).*(back_proj(ind)./(back_proj0(ind)+beta*image(ind)+beta*U(ind)));      % ������ֵ������ؽ�
    k=k+1;
    
    % ����betaֵ  ����Ӧ����betaֵ
    son = dot(delta,delta);
    mother = dot(U,U);
    beta = log(2*(son/(mother+0.1))+1);
    % ����ÿ�ε����Ľ��
    result_save = reshape(F,N,N)';
    writebin(['.\result\adaptive\adptive_' num2str(k) '.bin'],result_save);

    %% ���һ�ε�������һ���˲����� %%
    %     P = 9;
    %     hs = 0.01;
    %     hr = 10;
    %     filt_F = BilaFilter( image,P,hs,hr );
    %     F = reshape(filt_F,N2,1);
    %     J= ordfilt2(result_save,5,ones(3,4));% ���ж�άͳ��˳�����
    %     F = reshape(J,N2,1);
%     h=fspecial('sobel') ;
    %hsize=10;
    %sigma=5;
    %h=fspecial('gaussian',hsize,sigma);
%     h=fspecial('laplacian');
    %outImg = imfilter(result_save,h)';
    %F = reshape(outImg,N2,1);    
end
end

