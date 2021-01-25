function E = Func_Fit(Images,U,r)
%% definition
% -Images, MxNxC imatrix
% -U, MxNxb eigen matrix
% -r, r of Omega
% -E, MxN matrix, root residual
%% prepare
R = 2*r+1;
[M,N,C] = size(Images);
Num = M*N;

E = zeros(M,N,C);
% 1 eigen
U = cat(3,U,ones(M,N));
[~,~,B] = size(U);
BNum = B*Num;
%% prepare for fast inverse iteration
i_sp = repmat(reshape((1:BNum)',[B,Num]),[B,1]);
j_sp = repmat((1:BNum),[B,1]);
i_sp = i_sp(:);
j_sp = j_sp(:);
%% fit every page
for iter = 1 : C
    % prepare II
    I = Images(:,:,iter);
    mII = imboxfilt(I.*I,R); 
    mIIv = reshape(mII,[1,1,Num]);
    % prepare UU
    Uv  = reshape(permute(U,[3,1,2]),[B,1,Num]);       %  Uv = B,1,Num
    UU  = bsxfun(@times,Uv,permute(Uv,[2,1,3])); 
    UU  = reshape(permute(UU,[3,1,2]),[M,N,B*B]);      %  UU = M,N,BxB
    mUU = imboxfilt(UU,R);
    sUU = R*R*reshape(permute(mUU,[3,1,2]),[B,B,Num]); % mUU = B,B,Num
    % prepare UI
    UI = U.*repmat(I,[1,1,B]);  
    mUI = imboxfilt(UI,R);
    sUIv = R*R*reshape(permute(mUI,[3,1,2]),[B,1,Num]);
    % inverse, dense to sparse, sparse to dense
    C = sUU+0.00000001*repmat(eye(B),[1,1,Num]);
    C_sp = sparse(i_sp,j_sp,C(:),BNum,BNum); 
    d_sp = sparse((1:BNum)',ones(BNum,1),sUIv(:),BNum,1);
    a_sp = full(C_sp\d_sp);
    a = reshape(a_sp,[B,1,Num]);
    % calc residual
    e = sum(sum(repmat(a,[1,B,1]).*sUU,1).*permute(a,[2,1,3]),2)/(R*R)+...
        mIIv-2*sum(a.*sUIv,1)/(R*R);
    e(e<0) = 0;
    E(:,:,iter) = reshape(e,[M,N]);
end
E = E.^0.5;
E = mean(E,3);
