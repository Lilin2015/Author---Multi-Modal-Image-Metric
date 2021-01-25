function D = Func_DAISY(Img)
    %% prepare
    if size(Img,3)~=1, Img = rgb2gray(Img); end
    [m,n] = size(Img);
    sigma = 15/6*[1,2,3,4];
    %% calc G
    DX = [diff(Img,1,2),zeros(m,1)];
    DY = [diff(Img,1,1);zeros(1,n)];
    G(:,:,1) = -DY; G(:,:,2) = (-DY+DX)/1.414;
    G(:,:,3) =  DX; G(:,:,4) = ( DY+DX)/1.414;
    G(:,:,5) =  DY; G(:,:,6) = ( DY-DX)/1.414;
    G(:,:,7) = -DX; G(:,:,8) = (-DY-DX)/1.414;
    G(G<0) = 0;
    %% calc H
    H0 = imgaussfilt(G,sigma(1)); H0 = normalize(H0,3);
    H1 = imgaussfilt(G,sigma(2)); H1 = normalize(H1,3);
    H2 = imgaussfilt(G,sigma(3)); H2 = normalize(H2,3);
    H3 = imgaussfilt(G,sigma(4)); H3 = normalize(H3,3);
    H(:,:,:,1) = H0;
    H(:,:,:,2) = H1;
    H(:,:,:,3) = H2;
    H(:,:,:,4) = H3;
    %% calc D
    H(isnan(H)) = 0;
    D(:,:,:,1) = H(:,:,:,1);
    for k = 0 : 2
        d = (k+1)*5;
        D(:,:,:,8*k+2) = imtranslate(H(:,:,:,k+2),[ 0  d],'OutputView','same','FillValues',0);
        D(:,:,:,8*k+3) = imtranslate(H(:,:,:,k+2),[ d  0],'OutputView','same','FillValues',0);
        D(:,:,:,8*k+4) = imtranslate(H(:,:,:,k+2),[ 0 -d],'OutputView','same','FillValues',0);
        D(:,:,:,8*k+5) = imtranslate(H(:,:,:,k+2),[-d  0],'OutputView','same','FillValues',0);
        D(:,:,:,8*k+6) = imtranslate(H(:,:,:,k+2),[ d  d]/1.414,'OutputView','same','FillValues',0);
        D(:,:,:,8*k+7) = imtranslate(H(:,:,:,k+2),[-d  d]/1.414,'OutputView','same','FillValues',0);
        D(:,:,:,8*k+8) = imtranslate(H(:,:,:,k+2),[-d -d]/1.414,'OutputView','same','FillValues',0);
        D(:,:,:,8*k+9) = imtranslate(H(:,:,:,k+2),[ d -d]/1.414,'OutputView','same','FillValues',0);
    end
end

