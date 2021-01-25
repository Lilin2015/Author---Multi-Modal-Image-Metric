%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to produce "DASC spectral data.mat", "DASC stereo
% data.mat", "DASC stereo example.mat".
% We have not provide these offline data on Github, because they cost large
% memory (about 5 GB).
% If you want to reproduce our experimental results, you need to:
% 1. dowload DASC toolbox, "https://github.com/seungryong/DASC"
% 2. run this script
% 3. delete the empty cookies of DASC
% 4. run "metricComparison_spectral", "stereoExample_DASC"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set DASC configuration
M_half = 15;
N_half = 2;
epsil = 0.09;
downSize = 1;
sigma_s = 2;
sigma_r = 0.2;
iter = 1;

load rp1_Middlebury.mat;
load rp2_Middlebury.mat;
%% build "DASC spectral data.mat"
DASC = zeros(51,61,128,6,155);
DASCnoisy = DASC;
for i = 1 : 155
    i
    for j = 1 : 6
        G = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(j),'.png')));
        Gnoisy = im2double(imread(strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(i),'_',num2str(j),'_noisy.png')));
        DASC(:,:,:,j,i) = mexDASC_GF(G,M_half,N_half,rp1,rp2,epsil,downSize);
        DASCnoisy(:,:,:,j,i) = mexDASC_GF(Gnoisy,M_half,N_half,rp1,rp2,epsil,downSize);
    end
end
save('DASC spectral data','DASC','DASCnoisy');
%% build "DASC stereo data.mat"
DASC = zeros(51,61,128,2,115);
DASCnoisy = DASC;
for i = 1 : 115
    i
    for j = 1 : 2
        G = im2double(imread(strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(i),'_',num2str(j),'.png')));
        Gnoisy = im2double(imread(strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(i),'_',num2str(j),'_noisy.png')));
        DASC(:,:,:,j,i) = mexDASC_GF(G,M_half,N_half,rp1,rp2,epsil,downSize);
        DASCnoisy(:,:,:,j,i) = mexDASC_GF(Gnoisy,M_half,N_half,rp1,rp2,epsil,downSize);
    end
end
save('DASC stereo data','DASC','DASCnoisy');
%% build "DASC stereo example.mat"
 Limg = im2double(imread(strcat(pwd,'/Fig_stereo/left.jpg')));
 Rimg = im2double(imread(strcat(pwd,'/Fig_stereo/right.jpg')));
 DASC_L = mexDASC_GF(Limg,M_half,N_half,rp1,rp2,epsil,downSize);
 DASC_R = mexDASC_GF(Rimg,M_half,N_half,rp1,rp2,epsil,downSize);
 save('DASC stereo example','DASC_L','DASC_R');