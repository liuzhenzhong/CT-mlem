function F = medfuncAdaptMap( W_ind,W_dat,N,F0,P,irt_num ,beta)
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
        delta(ii) = P(ii) - proj;  % 投影与前向投影的差值向量
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
            back_proj0(index) = back_proj0(index)+data(ind).';%+ beta*U(index);%全1向量的反投影
        end
    end
    % =========图像更新一次============== %
    ind = back_proj0>0;
%     F(ind) = F(ind).*(back_proj(ind)./back_proj0(ind));   % 普通的ml_em算法
    F(ind) = F(ind).*(back_proj(ind)./(back_proj0(ind)+beta*U(ind)));    % 基于马尔科夫随机场先验知识的重建
%     F(ind) = F(ind).*(back_proj(ind)./(back_proj0(ind)+beta*image(ind)+beta*U(ind)));      % 基于中值先验的重建
    k=k+1;
    
    % 更新beta值  自适应更新beta值
    son = dot(delta,delta);
    mother = dot(U,U);
    beta = log(2*(son/(mother+0.1))+1);
    % 保存每次迭代的结果
    result_save = reshape(F,N,N)';
    writebin(['.\result\adaptive\adptive_' num2str(k) '.bin'],result_save);

    %% 完成一次迭代进行一次滤波处理 %%
    %     P = 9;
    %     hs = 0.01;
    %     hr = 10;
    %     filt_F = BilaFilter( image,P,hs,hr );
    %     F = reshape(filt_F,N2,1);
    %     J= ordfilt2(result_save,5,ones(3,4));% 进行二维统计顺序过滤
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

