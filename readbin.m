function  a = readbin(file)
fid = fopen(file,'rb');
fseek(fid,-20,1);
 fmin = fread(fid,1,'float');
 fmax = fread(fid,1,'float');
 width = fread(fid,1,'long');
 height = fread(fid,1,'long');
 depth = fread(fid,1,'long');
fclose(fid);

fid = fopen(file,'rb');
a = fread(fid,[width height],'float');%size [M��N]�������ݵ�M��N�ľ����У����ݰ��д�ţ� inf�������ļ�
a = a';
fclose(fid);