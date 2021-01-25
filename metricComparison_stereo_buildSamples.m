%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script shows how we build testing samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error('Exising samples will be covered! The results of DASC and L2-Net must be reformed! If you are sure, comment out this line.')

clear
close all

maxShift = 5;           % the maximal shift in experiment
R = 25;                 % the radius of \Omega
SideL = maxShift + R;   % it affects the size of cropped subImages

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% crop the original samples
Img = zeros(512,512,2);
G   = zeros(512,512,2);
for i = 1 : 23
    fprintf('\n%d/23',i); pause(0.0001);
    % read images
    temp = im2double(imread(...
            strcat(pwd,'/Fig_metricComparison/Stereo Origin/',num2str(i),...
            '_1.png')));
    if size(temp,3)~=1, temp=rgb2gray(temp); end
    Img(:,:,1) = imresize(temp,[512 512]);
    temp = pfmread(...
            strcat(pwd,'/Fig_metricComparison/Stereo Origin/',num2str(i),...
            '_2.pfm'));
    if size(temp,3)~=1, temp=rgb2gray(temp); end
    Img(:,:,2) = imresize(temp/max(temp(:)),[512 512]);
    % crop
    [M,N,~] = size(Img);
    minM = 1+SideL; minN = minM;    % rows range 
    maxM = M-SideL; maxN = N-SideL; % column range
    [G(:,:,1),~] = imgradient(Img(:,:,1));
    [G(:,:,2),~] = imgradient(Img(:,:,2));
    G = min(G,[],3);
    G(1:SideL,:) = 0; % remove padding
    G(M-SideL+1:M,:) = 0;
    G(:,1:SideL) = 0;
    G(:,N-SideL+1:N) = 0;
    % write
    for k = 1 : 5
        [m,n] = find(G==max(G(:))); % find the most sharp region
        m=m(1); n=n(1);
        ImgNoise = imnoise(Img,'gaussian',0,0.001);
        imwrite(Img(m-R:m+R,n-SideL:n+SideL,1),strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(5*(i-1)+k),'_1.png'));
        imwrite(Img(m-R:m+R,n-SideL:n+SideL,2),strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(5*(i-1)+k),'_2.png'));
        imwrite(ImgNoise(m-R:m+R,n-SideL:n+SideL,1),strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(5*(i-1)+k),'_1_noisy.png'));
        imwrite(ImgNoise(m-R:m+R,n-SideL:n+SideL,2),strcat(pwd,'/Fig_metricComparison/Stereo Samples/',num2str(5*(i-1)+k),'_2_noisy.png'));
        G(m-R:m+R,n-SideL:n+SideL,:)=0;
    end
end
fprintf('\n');