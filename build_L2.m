%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to produce "L2Net spectral data.mat", "L2Net stereo
% data.mat", "L2Net stereo example.mat".
% We have not provide these offline data on Github, because they cost large
% memory (more than 5 GB).
% If you want to reproduce our experimental results, you need to 
% (this will take a very long time, about 8 hours, if you can not run matcovnet by GPU):
% 1. dowload L2Net toolbox, "https://github.com/yuruntian/L2-Net"
% 2. install matconvnet
% 3. delete the empty cookies of L2Net
% 4. run "metricComparison_spectral", "stereoExample_DASC"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set L2Net configuration
run(fullfile('vl_setupnn.m')) ;
rootPath = pwd;
trainSet = 'LIB';%'YOS','ND','HP'
flagAug = 1;
flagGPU = 1;
batchSize = 1000;
%% build "L2Net spectral data.mat"
L2 = zeros(4,256,7,6,155);
L2noisy = zeros(4,256,7,6,155);
for i = 1 : 155
    i
    for j = 1 : 6
        Img  = imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(j),'.png'));
        ImgN = imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(j),'_noisy.png'));
        iter = 0;
        for s = [-5,-2,-1,0,1,2,5]
            iter = iter + 1;
            subImg = imresize(Img(:,6+s:56+s),[64,64]);
            subImgN = imresize(ImgN(:,6+s:56+s),[64,64]);
            L2(:,:,iter,j,i) = cal_L2Net_des(rootPath,trainSet,flagCS,flagAug,single(subImg),batchSize,flagGPU);
            L2noisy(:,:,iter,j,i) = cal_L2Net_des(rootPath,trainSet,flagCS,flagAug,single(subImgN),batchSize,flagGPU);
        end
    end
end
save('L2Net spectral data','L2','L2noisy');
%% build "L2Net stereo data.mat"
L2 = zeros(4,256,7,2,115);
L2noisy = zeros(4,256,7,2,115);
for i = 1 : 115
    i
    for j = 1 : 2
        Img  = imread(strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(i),'_',num2str(j),'.png'));
        ImgN = imread(strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(i),'_',num2str(j),'_noisy.png'));
        iter = 0;
        for s = [-5,-2,-1,0,1,2,5]
            iter = iter + 1;
            subImg = imresize(Img(:,6+s:56+s),[64,64]);
            subImgN = imresize(ImgN(:,6+s:56+s),[64,64]);
            L2(:,:,iter,j,i) = cal_L2Net_des(rootPath,trainSet,flagCS,flagAug,single(subImg),batchSize,flagGPU);
            L2noisy(:,:,iter,j,i) = cal_L2Net_des(rootPath,trainSet,flagCS,flagAug,single(subImgN),batchSize,flagGPU);
        end 
    end
end
save('L2Net stereo data','L2','L2noisy');
%% build "L2Net stereo example.mat"
% it takes me 7 hours to build a L2-Net dense descriptor for this sample
% if you do not have 7 hours, you can terminate this process, and copy the last sentence to save
% current results, and continue next time.
L = imread(strcat(pwd,'/Fig_stereo/left.jpg'));
R = imread(strcat(pwd,'/Fig_stereo/right.jpg'));
L2_L = zeros(4,256,375,450);
L2_R = zeros(4,256,375,450);
load('L2Net stereo example.mat');
for im = 16 : 375-15
    for in = 16 : 450-15
        [im,in]
        if any(L2_L(:,:,im,in))
            continue;
        end
        subL = imresize(L(im-15:im+15,in-15:in+15),[64 64]);
        subR = imresize(R(im-15:im+15,in-15:in+15),[64 64]);
        L2_L(:,:,im,in) = cal_L2Net_des(rootPath,trainSet,flagCS,flagAug,single(subL),batchSize,flagGPU);
        L2_R(:,:,im,in) = cal_L2Net_des(rootPath,trainSet,flagCS,flagAug,single(subR),batchSize,flagGPU);
    end
end
save('L2Net stereo example','L2_L','L2_R');