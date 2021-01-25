function [MaskL,MaskR] = Func_crossCheck(DispL,DispR)
%% ������֤, ����stereo images��Ч
% DispL(x,y)=k, ��ʾ��ͼ(x,y)��������ͼ������Ϊ(x+k,y)
% DispR(x,y)=k, ��ʾ��ͼ(x,y)��������ͼ������Ϊ(x-k,y)
% MaskL, MaskR, 1��ʾ������֤ͨ��, ����Ϊ0

[M,N] = size(DispL);

Midx = repmat((1:M)',[1,N]);
Nidx = repmat((1:N),[M,1]);
%% left
Nidx_toL = max(1,Nidx-DispL);
idx_map = sub2ind([M,N],Midx(:),Nidx_toL(:));
MaskL = reshape(abs(DispR(idx_map)-DispL(:))<1,[M,N]);
%% right
Nidx_toR = min(N,Nidx+DispR);
idx_map = sub2ind([M,N],Midx(:),Nidx_toR(:));
MaskR = reshape(abs(DispL(idx_map)-DispR(:))<1,[M,N]);

end

%m=90
%n=85