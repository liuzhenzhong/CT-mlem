function F = medfuncMlem( W_ind,W_dat,N,F0,P,irt_num )
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
    R = zeros(M,1);
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
            back_proj0(index) = back_proj0(index)+data(ind).';%ȫ1�����ķ�ͶӰ
        end
    end
    % =========ͼ�����һ��============== %
    ind = back_proj0>0;
    F(ind) = F(ind).*(back_proj(ind)./back_proj0(ind));
    k=k+1;
    
    % ����ÿ�ε����Ľ��
    result_save = reshape(F,N,N)';
    writebin(['.\result\mlem_' num2str(k) '.bin'],result_save); 
end
end

