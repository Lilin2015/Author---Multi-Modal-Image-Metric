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
Img = zeros(512,512,6);
G   = zeros(512,512,6);
for i = 1 : 31
    fprintf('\n%d/31',i); pause(0.0001);
    % read original samples
    for k = 1 : 6
        temp = im2double(imread(...
            strcat(pwd,'/Fig_metricComparison/Multi-Spectral Origin/',num2str(i),...
            '_',num2str(k),'.png')));
        if size(temp,3)~=1, temp=rgb2gray(temp); end
        Img(:,:,k) = temp;
    end
    % crop
    [M,N,~] = size(Img);
    minM = 1+SideL; minN = minM;    % rows range 
    maxM = M-SideL; maxN = N-SideL; % column range
    for k = 1 : 6
        [G(:,:,k),~] = imgradient(Img(:,:,k));
    end
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
        for t = 1 : 6
            subImg = Img(m-R:m+R,n-SideL:n+SideL,t);
            imwrite(subImg,strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(5*(i-1)+k),'_',num2str(t),'.png'));
            subImgNoise = ImgNoise(m-R:m+R,n-SideL:n+SideL,t);
            imwrite(subImgNoise,strcat(pwd,'/Fig_metricComparison/Multi-Spectral Samples/',num2str(5*(i-1)+k),'_',num2str(t),'_noisy.png'));
        end
        G(m-R:m+R,n-SideL:n+SideL,:)=0;
    end
end
fprintf('\n');