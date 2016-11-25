function F = medfuncMlem( W_ind,W_dat,N,F0,P,irt_num )
%UNTITLED5 此处显示有关此函数的摘要
%   输入参数：
%   输出参数：
N2 = N^2;
M = length(P); % 所有射线条数
F = F0;
k=0;
while(k<irt_num)
    disp(['共' num2str(irt_num) '次迭代，第' num2str(k+1) '次迭代']);
    %% ==求实际投影值和估计投影值的比值向量=== %%
    R = zeros(M,1);
    for ii=1:M
        % 如果射线不通过任何像素，不做计算
        if any(W_ind(ii,:)) == 0
            continue;
        end
        w = zeros(1,N2);%系统矩阵的一个行向量
        for jj=1:2*N
            m=W_ind(ii,jj);
            if m>0&&m<=N2
                w(m)=W_dat(ii,jj);
            end
        end
        proj = w*F;%前向投影
        R(ii)=P(ii)./proj;%比值向量，一条射线一个比值
    end
    % ====== 求比值向量与全1向量的反投影 ========= %
    back_proj = zeros(N2,1);back_proj0=zeros(N2,1);
    for ii=1:M %以行的遍历来进行反投影计算
        label = W_ind(ii,:);
        data = W_dat(ii,:);
        if any(label)~=0
            ind = label>0;
            index = label(ind);
            back_proj(index) = back_proj(index)+R(ii)*data(ind).';%比值向量的反投影
            back_proj0(index) = back_proj0(index)+data(ind).';%全1向量的反投影
        end
    end
    % =========图像更新一次============== %
    ind = back_proj0>0;
    F(ind) = F(ind).*(back_proj(ind)./back_proj0(ind));
    k=k+1;
    
    % 保存每次迭代的结果
    result_save = reshape(F,N,N)';
    writebin(['.\result\mlem_' num2str(k) '.bin'],result_save); 
end
end

