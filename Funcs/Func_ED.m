function score = Func_ED(A,B)
    edgeA = edge(A,'Canny');
    edgeB = edge(B,'Canny');
    [~,GA] = imgradient(A);
    [~,GB] = imgradient(B);
    
    Ggap = abs(GA-GB);
    Ggap(Ggap>180)=360-Ggap(Ggap>180);
    Gsame = (Ggap<(360/16));
    
    matchPts = and(and(Gsame,edgeA),edgeB);
    if matchPts==0
        score = 0;
    else
        score = sum(matchPts(:))/(sum(edgeB(:)).^0.5);
    end
end

