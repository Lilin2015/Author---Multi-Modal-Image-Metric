%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script shows how we build the stereo results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,strcat(pwd,'/Funcs'))

imgL = im2double(imread(strcat(pwd,'/Fig_stereo/left.jpg')));
imgR = im2double(imread(strcat(pwd,'/Fig_stereo/right.jpg')));
load('L2Net stereo example.mat');

feaL = permute(L2_L,[3,4,1,2]);
feaR = permute(L2_R,[3,4,1,2]);
feaL(isnan(feaL)) = 0;
feaR(isnan(feaR)) = 0;

CL = zeros(375,450,60);
CR = zeros(375,450,60);

for i = 0 : 60
    i
    feaL_Mov = imtranslate(feaL,[-i,0],'OutputView','same','FillValues',0);
    feaR_Mov = imtranslate(feaR,[i,0],'OutputView','same','FillValues',0);
    
    E_L = mean(vecnorm(feaR_Mov-feaL,2,3),4);
    E_R = mean(vecnorm(feaL_Mov-feaR,2,3),4);
    
    smoothE_L = imguidedfilter(E_L,imgL,'NeighborhoodSize',[4 4],'DegreeOfSmoothing',0.0001);
    smoothE_R = imguidedfilter(E_R,imgR,'NeighborhoodSize',[4 4],'DegreeOfSmoothing',0.0001);
    
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
imwrite(ind2rgb(gray2ind(LM./max(LM(:)),255),jet(255)),strcat(pwd,'/Fig_stereo/stereo_L2.png'));