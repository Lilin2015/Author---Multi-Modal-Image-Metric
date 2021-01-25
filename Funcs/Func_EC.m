function E = Func_EC(A,B,K,r)        
    Imgs = cat(3,A,B); 
    for i = 1 : size(Imgs,3)
        Imgs(:,:,i) = Imgs(:,:,i) - mean(mean(Imgs(:,:,i)));
    end
    [M,N,C] = size(Imgs);
    X = reshape(Imgs,[M*N,C]);
    [U,~,~] = svds(X,K); 
    Base = reshape(U,[M,N,K]);
    E = Func_Fit(Imgs,Base,r);
end

