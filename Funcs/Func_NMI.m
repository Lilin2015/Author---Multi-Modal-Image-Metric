% calculate normalized mutual information between image A and B
function MI = Func_NMI( A, B )
    %% Pxy
    Xedges = (0:1:255)/255;
    Yedges = (0:1:255)/255;
    [H,~,~] = histcounts2(A,B,Xedges,Yedges);
    Pxy = H/sum(H(:));
    %% Px and Py
    Px = sum(Pxy,2);
    Py = sum(Pxy,1);
    %% entropy
    Hxy = Pxy.*log(Pxy./(Px.*Py)); 
    Hx  = -Px.*log(Px);
    Hy  = -Py.*log(Py);
    Hxy = sum(Hxy(:),'omitnan');
    Hx = sum(Hx(:),'omitnan');
    Hy = sum(Hy(:),'omitnan');
    %% Normalized Mutual information
    MI = 2*Hxy/(Hx+Hy);
end