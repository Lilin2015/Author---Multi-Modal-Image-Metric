%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function test metrics based on the samples
% the results are in workspace, named "CorrectRatio_XX"
% the cookies are in "Fig_metricComparison/XXX/XXX Cookies"
% if you want to update the results, delete the cookies
% if you want to test a specific method, delete its cookies only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
path(path,strcat(pwd,'/Funcs'))

% set parameters
pair  = [3,4];              % 1-6 means 400,460,520,580,640,700nm repectively
Noise = 0;                  % 0- test on noise-free samples, 1- test on noisy samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify the following parameters, unless you have read though
% all the scripts, and know what a chain reaction it might trigger
Shift = [-5,-2,-1,0,1,2,5];
Patch_r = 25;              
Num = 155;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Record_MI  = zeros(Num,size(Shift,2));
Record_CC  = Record_MI;
Record_DAISY = Record_MI;
Record_DASC = Record_MI;
Record_ED = Record_MI;
Record_L2 = Record_MI;
Record_EC = Record_MI;
Record_MI6  = Record_MI;
Record_CC6  = Record_MI;
Record_DAISY6 = Record_MI;
Record_DASC6 = Record_MI;
Record_ED6 = Record_MI;
Record_L26 = Record_MI;
Record_EC6 = Record_MI;

noiseStr = '';
if Noise == 1
    noiseStr = '_noisy';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test MI
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_MI.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test MI %d/',num2str(Num)),i); pause(0.0001);
        % read images
        ImgRef = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(2)),noiseStr,'.png')));
        ImgMov = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(1)),noiseStr,'.png')));
        % shift the images for testing
        center = [(size(ImgRef,1)+1)/2,(size(ImgRef,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = ImgMov(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s);
            ImgRefd = ImgRef(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r);
            
            Record_MI(i,iter_s)  = Func_NMI(ImgRefd,ImgMovd);
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_MI.mat'),'Record_MI');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test MI6
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_MI6.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test MI6 %d/',num2str(Num)),i); pause(0.0001);
        % read images
        Img(:,:,1) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_1',noiseStr,'.png')));
        Img(:,:,2) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_2',noiseStr,'.png')));
        Img(:,:,3) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_3',noiseStr,'.png')));
        Img(:,:,4) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_4',noiseStr,'.png')));
        Img(:,:,5) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_5',noiseStr,'.png')));
        Img(:,:,6) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_6',noiseStr,'.png')));
        % shift the images for testing
        center = [(size(Img,1)+1)/2,(size(Img,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s,1:3);
            ImgRefd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r,4:6);

            for t = 1 : 3
                MItemp(t) = Func_NMI(ImgRefd(:,:,t),ImgMovd(:,:,t));
            end
            
            Record_MI6(i,iter_s) = mean(MItemp);
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_MI6.mat'),'Record_MI6');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test CC
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_CC.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test CC %d/',num2str(Num)),i); pause(0.0001);
        % read image pair (left is fixed, right is floating)
        ImgRef = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(2)),noiseStr,'.png')));
        ImgMov = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(1)),noiseStr,'.png')));
        % shift the images for testing
        center = [(size(ImgRef,1)+1)/2,(size(ImgRef,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = ImgMov(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s);
            ImgRefd = ImgRef(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r);
            
            Record_CC(i,iter_s)  = corr2(ImgRefd,ImgMovd);
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_CC.mat'),'Record_CC');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test CC6
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_CC6.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test CC6 %d/',num2str(Num)),i); pause(0.0001);
        % read images
        Img(:,:,1) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_1',noiseStr,'.png')));
        Img(:,:,2) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_2',noiseStr,'.png')));
        Img(:,:,3) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_3',noiseStr,'.png')));
        Img(:,:,4) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_4',noiseStr,'.png')));
        Img(:,:,5) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_5',noiseStr,'.png')));
        Img(:,:,6) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_6',noiseStr,'.png')));
        % shift the images for testing
        center = [(size(Img,1)+1)/2,(size(Img,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s,1:3);
            ImgRefd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r,4:6);
            
            for t = 1:3
                CCtemp(t) = corr2(ImgRefd(:,:,t),1-ImgMovd(:,:,t));
            end
            
            Record_CC6(i,iter_s) = mean(CCtemp);
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_CC6.mat'),'Record_CC6');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test DAISY
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DAISY.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test DAISY %d/',num2str(Num)),i); pause(0.0001);
        % read image pair (left is fixed, right is floating)
        ImgRef = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(2)),noiseStr,'.png')));
        ImgMov = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(1)),noiseStr,'.png')));
        % shift the images for testing
        center = [(size(ImgRef,1)+1)/2,(size(ImgRef,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = ImgMov(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s);
            ImgRefd = ImgRef(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r);
            
            DAISY_Mov = Func_DAISY(ImgMovd);
            DAISY_Ref = Func_DAISY(ImgRefd);
            Record_DAISY(i,iter_s)  = mean(vecnorm(DAISY_Mov(26,26,:,:)-DAISY_Ref(26,26,:,:),2,3));
        end 
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DAISY.mat'),'Record_DAISY');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test DAISY6
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DAISY6.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test DAISY6 %d/',num2str(Num)),i); pause(0.0001);
        % read images
        Img(:,:,1) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_1',noiseStr,'.png')));
        Img(:,:,2) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_2',noiseStr,'.png')));
        Img(:,:,3) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_3',noiseStr,'.png')));
        Img(:,:,4) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_4',noiseStr,'.png')));
        Img(:,:,5) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_5',noiseStr,'.png')));
        Img(:,:,6) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_6',noiseStr,'.png')));
        % shift the images for testing
        center = [(size(Img,1)+1)/2,(size(Img,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s,1:3);
            ImgRefd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r,4:6);
            
            for t = 1:3
                DAISY_Mov = Func_DAISY(ImgMovd(:,:,t));
                DAISY_Ref = Func_DAISY(ImgRefd(:,:,t));
                DAISYtemp(t) = mean(vecnorm(DAISY_Mov(26,26,:,:)-DAISY_Ref(26,26,:,:),2,3));
            end
            Record_DAISY6(i,iter_s) = mean(DAISYtemp);
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DAISY6.mat'),'Record_DAISY6');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test DASC
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DASC.mat'),'file')==0
    load('DASC spectral data.mat');
    for i = 1 : Num
        fprintf(strcat('\n test DASC %d/',num2str(Num)),i); pause(0.0001);
        % read DASC_data
        if Noise == 0
            DASC_Ref = DASC(:,:,:,pair(2),i);
            DASC_Mov = DASC(:,:,:,pair(1),i);
        else
            DASC_Ref = DASCnoisy(:,:,:,pair(2),i);
            DASC_Mov = DASCnoisy(:,:,:,pair(1),i);
        end
        % shift the features for testing
        center = [(size(DASC_Ref,1)+1)/2,(size(DASC_Ref,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            DASC_Movd = DASC_Mov(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s,:);
            DASC_Refd = DASC_Ref(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r,:);
            Record_DASC(i,iter_s) = mean(mean(vecnorm(DASC_Movd(17:35,17:35,:)-DASC_Refd(17:35,17:35,:),2,3),'omitnan'),'omitnan');
        end 
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DASC.mat'),'Record_DASC');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test DASC6
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DASC6.mat'),'file')==0
    load('DASC spectral data.mat');
    for i = 1 : Num
        fprintf(strcat('\n test DASC6 %d/',num2str(Num)),i); pause(0.0001);
        % read DASC_data
        if Noise == 0
            DASC_Mov = DASC(:,:,:,1:3,i);
            DASC_Ref = DASC(:,:,:,4:6,i);
        else
            DASC_Mov = DASCnoisy(:,:,:,1:3,i);
            DASC_Ref = DASCnoisy(:,:,:,4:6,i);
        end
        % shift the features for testing
        center = [(size(DASC_Ref,1)+1)/2,(size(DASC_Ref,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            for t = 1:3
                DASC_Movd = DASC_Mov(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s,:,t);
                DASC_Refd = DASC_Ref(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r,:,t);
                temp_DASC(t)  = mean(mean(vecnorm(DASC_Movd(17:35,17:35,:)-DASC_Refd(17:35,17:35,:),2,3),'omitnan'),'omitnan');
            end
            Record_DASC6(i,iter_s) = mean(temp_DASC);
        end 
        
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DASC6.mat'),'Record_DASC6');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test ED
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_ED.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test ED %d/',num2str(Num)),i); pause(0.0001);
        % read image pair (left is fixed, right is floating)
        ImgRef = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(2)),noiseStr,'.png')));
        ImgMov = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(1)),noiseStr,'.png')));
        % shift the images for testing
        center = [(size(ImgRef,1)+1)/2,(size(ImgRef,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = ImgMov(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s);
            ImgRefd = ImgRef(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r);
            
            Record_ED(i,iter_s)  = Func_ED(ImgRefd,ImgMovd);
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_ED.mat'),'Record_ED');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test ED6
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_ED6.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test ED6 %d/',num2str(Num)),i); pause(0.0001);
        % read images
        Img(:,:,1) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_1',noiseStr,'.png')));
        Img(:,:,2) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_2',noiseStr,'.png')));
        Img(:,:,3) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_3',noiseStr,'.png')));
        Img(:,:,4) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_4',noiseStr,'.png')));
        Img(:,:,5) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_5',noiseStr,'.png')));
        Img(:,:,6) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_6',noiseStr,'.png')));
        % shift the images for testing
        center = [(size(Img,1)+1)/2,(size(Img,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s,[1,2,3]);
            ImgRefd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r,4:6);
            
            for t = 3
                ED6temp(t) = Func_ED(ImgRefd(:,:,t),ImgMovd(:,:,t));
            end
            
            Record_ED6(i,iter_s) = mean(ED6temp);
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_ED6.mat'),'Record_ED6');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test L2
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_L2.mat'),'file')==0
    load('L2Net spectral data');
    for i = 1 : Num
        fprintf(strcat('\n test 2L %d/',num2str(Num)),i); pause(0.0001);
        % read feature of image pair
        for t = 1:7
            fprintf('.');
            if Noise == 0
                L2Movd = L2(:,:,t,pair(1),i);
                L2Refd = L2(:,:,t,pair(2),i);
            else
                L2Movd = L2noisy(:,:,t,pair(1),i);
                L2Refd = L2noisy(:,:,t,pair(2),i);
            end
            Record_L2(i,t) = mean(vecnorm(L2Movd-L2Refd,2,1));
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_L2.mat'),'Record_L2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test L26
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_L26.mat'),'file')==0
    load('L2Net spectral data');
    for i = 1 : Num
        fprintf(strcat('\n test 2L %d/',num2str(Num)),i); pause(0.0001);
        % read feature of images
        for t = 1:7
            fprintf('.');
            if Noise == 0
                L2Movd = L2(:,:,t,1:3,i);
                L2Refd = L2(:,:,t,4:6,i);
            else
                L2Movd = L2noisy(:,:,t,1:3,i);
                L2Refd = L2noisy(:,:,t,4:6,i);
            end
            Record_L26(i,t) = mean(mean(vecnorm(L2Movd-L2Refd,2,1)));
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_L26.mat'),'Record_L26');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test EC
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_EC.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test EC %d/',num2str(Num)),i); pause(0.0001);
        % read image pair
        ImgRef = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(2)),noiseStr,'.png')));
        ImgMov = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(pair(1)),noiseStr,'.png')));
        % shift the images for testing
        center = [(size(ImgRef,1)+1)/2,(size(ImgRef,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = ImgMov(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s);
            ImgRefd = ImgRef(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r);
            
            E = Func_EC(ImgRefd,ImgMovd,1,3);
            Record_EC(i,iter_s) = mean(E(:));
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_EC.mat'),'Record_EC');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test EC-6 channel
if exist(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_EC6.mat'),'file')==0
    for i = 1 : Num
        fprintf(strcat('\n test EC-6 %d/',num2str(Num)),i); pause(0.0001);
        % read images
        Img(:,:,1) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_1',noiseStr,'.png')));
        Img(:,:,2) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_2',noiseStr,'.png')));
        Img(:,:,3) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_3',noiseStr,'.png')));
        Img(:,:,4) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_4',noiseStr,'.png')));
        Img(:,:,5) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_5',noiseStr,'.png')));
        Img(:,:,6) = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_6',noiseStr,'.png')));
        % shift the images for testing
        center = [(size(Img,1)+1)/2,(size(Img,2)+1)/2];
        iter_s = 0;
        for s = Shift
            fprintf('.');
            iter_s = iter_s + 1;
            ImgMovd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r+s:center(2)+Patch_r+s,1:3);
            ImgRefd = Img(center(1)-Patch_r:center(1)+Patch_r,center(2)-Patch_r:center(2)+Patch_r,4:6);
            
            E = Func_EC(ImgRefd,ImgMovd,3,12);
            Record_EC6(i,iter_s) = mean(E(:));
        end
    end
    fprintf('\n');
    save(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_EC6.mat'),'Record_EC6');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statistic
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_MI.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_CC.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DAISY.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DASC.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_ED.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_L2.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_EC.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_MI6.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_CC6.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DAISY6.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_DASC6.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_ED6.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_L26.mat'));
load(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Cookies/Record_EC6.mat'));

[~,R_MI_opt] = max(Record_MI,[],2);
[~,R_CC_opt] = max(abs(Record_CC),[],2);
[~,R_DAISY_opt] = min(Record_DAISY,[],2);
[~,R_DASC_opt] = min(Record_DASC,[],2);
[~,R_ED_opt] = max(Record_ED,[],2);
[~,R_L2_opt] = max(Record_L2,[],2);
[~,R_EC_opt] = min(Record_EC,[],2);
[~,R_MI6_opt] = max(Record_MI6,[],2);
[~,R_CC6_opt] = max(abs(Record_CC6),[],2);
[~,R_DAISY6_opt] = min(Record_DAISY6,[],2);
[~,R_DASC6_opt] = min(Record_DASC6,[],2);
[~,R_ED6_opt] = max(Record_ED6,[],2);
[~,R_L26_opt] = max(Record_L26,[],2);
[~,R_EC6_opt] = min(Record_EC6,[],2);

correct = (size(Shift,2)+1)/2;
CorrectRatio_MI = sum(abs(R_MI_opt-correct)<1)/Num;
CorrectRatio_CC = sum(abs(R_CC_opt-correct)<1)/Num;
CorrectRatio_DAISY = sum(abs(R_DAISY_opt-correct)<1)/Num;
CorrectRatio_DASC = sum(abs(R_DASC_opt-correct)<1)/Num;
CorrectRatio_ED = sum(abs(R_ED_opt-correct)<1)/Num;
CorrectRatio_L2 = sum(abs(R_L2_opt-correct)<1)/Num;
CorrectRatio_EC = sum(abs(R_EC_opt-correct)<1)/Num;
CorrectRatio_MI6 = sum(abs(R_MI6_opt-correct)<1)/Num;
CorrectRatio_CC6 = sum(abs(R_CC6_opt-correct)<1)/Num;
CorrectRatio_DAISY6 = sum(abs(R_DAISY6_opt-correct)<1)/Num;
CorrectRatio_DASC6 = sum(abs(R_DASC6_opt-correct)<1)/Num;
CorrectRatio_ED6 = sum(abs(R_ED6_opt-correct)<1)/Num;
CorrectRatio_L26 = sum(abs(R_L26_opt-correct)<1)/Num;
CorrectRatio_EC6 = sum(abs(R_EC6_opt-correct)<1)/Num;