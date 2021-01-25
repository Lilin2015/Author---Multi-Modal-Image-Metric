%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script shows how we build the stereo results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,strcat(pwd,'/Funcs'))

imgL = im2double(imread(strcat(pwd,'/Fig_stereo/left.jpg')));
imgR = im2double(imread(strcat(pwd,'/Fig_stereo/right.jpg')));

CL = zeros(375,450,60);
CR = zeros(375,450,60);

for i = 0 : 60
    i
    imgL_Mov = imtranslate(imgL,[-i,0],'OutputView','same','FillValues',0);
    imgR_Mov = imtranslate(imgR,[i,0],'OutputView','same','FillValues',0);
    
    E_L = Func_EC(imgL,imgR_Mov,1,3);
    E_R = Func_EC(imgL_Mov,imgR,1,3);
    
    smoothE_L = imguidedfilter(E_L,imgL,'NeighborhoodSize',[9 9],'DegreeOfSmoothing',0.0001);
    smoothE_R = imguidedfilter(E_R,imgR,'NeighborhoodSize',[9 9],'DegreeOfSmoothing',0.0001);
    
    CL(:,:,i+1) = smoothE_L;
    CR(:,:,i+1) = smoothE_R;
end

[~,RoughDL]=min(CL,[],3);
[~,RoughDR]=min(CR,[],3);
[MaskL,MaskR] = Func_crossCheck(RoughDL-1,RoughDR-1);

LM = RoughDL-1;
LM(MaskL==0)=NaN;

GT = double(imread(strcat(pwd,'/Fig_stereo/GT.png')))/4;
ratio = sum(MaskL(:))/(375*450)                 % valid ratio
accuracy = sum(sum(abs(GT-LM)<=1))/sum(MaskL(:))   % accuracy

imshow(ind2rgb(gray2ind(LM./max(LM(:)),255),jet(255)));
imwrite(ind2rgb(gray2ind(LM./max(LM(:)),255),jet(255)),strcat(pwd,'/Fig_stereo/stereo_EC.png'));