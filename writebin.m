function writebin(file,a)
[height width] = size(a);
fid = fopen(file,'wb');
fwrite(fid,a','float');
k = zeros(1,123);
k = double(k);
fwrite(fid,k,'float');
k = min(a(:));
fwrite(fid,k,'float');
k = max(a(:));
fwrite(fid,k,'float');
fwrite(fid,width,'int');
fwrite(fid,height,'int');
Depth = 1;
fwrite(fid,Depth,'int');
fclose(fid);