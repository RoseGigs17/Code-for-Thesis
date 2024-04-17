% try with the 16 and 36 cm diameters first 

% % 16 cm % % 
% 30%
Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files16_30 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010003/D00030018","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0018", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020018"];

% 60%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files16_60 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010004/D00040018", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0018", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020018"];

% 90%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files16_90 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010005/D00050018", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0018", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020018"];

% % 21 cm % % 
% 30%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files21_30 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010003/D00030043","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0043", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020043"];

% 60%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files21_60 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010004/D00040043", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0043", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020043"];

% 90%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files21_90 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010005/D00050043", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0043", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020043"];

% % 26 cm % % 
% 30%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files26_30 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010003/D00030069","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0069", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020069"];

% 60%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files26_60 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010004/D00040069", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0069", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020069"];

% 90%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files26_90 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010005/D00050069", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0069", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020069"];

% % 31 cm % % 
% 30%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files31_30 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010003/D00030080","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0080", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020080"];

% 60%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files31_60 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010004/D00040080", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0080", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020080"];

% 90%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files31_90 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010005/D00050080", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0080", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020080"];


% % 36 cm % % 
% 30%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files36_30 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010003/D00030109","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0109", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020109"];

% 60%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files36_60 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010004/D00040109", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0109", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020109"];

% 90%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
Files36_90 = ["Mercury4/Neusoft AI Mercury4 Image Data/00010005/D00050109", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0109", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020109"];

%%%% attempt at making one giant data%%%%
diameter = [16,21,26,31,36];
Files30 = vertcat(Files16_30, Files21_30, Files26_30, Files31_30, Files36_30);
Files60 = vertcat(Files16_60, Files21_60, Files26_60, Files31_60, Files36_60);
Files90 = vertcat(Files16_90, Files21_90, Files26_90, Files31_90, Files36_90);

%Files = {Files30, Files60, Files90};

%large ROI
% ROI = [200,210,320,330]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
% ROI_display = [200, 210, 120, 120] %[xmin ymin xmax ymax]

% large roi no middle rod
% ROI = [190,280,330,330]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
% ROI_display = [190, 280, 140, 50]; %[xmin ymin xmax ymax]

%bone ROI
% ROI = [310,240,340,250];
% ROI_disp = [310,240,30,10];

%polystyrene ROI
% ROI = [300,305,330,315];
% ROI_disp = [300,305,30,10];

%iodine ROI
% ROI = [230,315,260,325];
% ROI_disp = [230,315,30,10];

%air ROI
% ROI = [205,255,235,265];
% ROI_disp = [205,255,30,10];

%water ROI
% ROI = [250,207,280,217];
% ROI_disp = [250,207,30,10];

%small ROI no middle rod
ROI = [270,266,320,316]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
ROI_display = [270, 266, 50, 50] %[xmin ymin xmax ymax]

%small ROI middle rod
% ROI = [250,266,300,316]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
% ROI_display = [250, 266, 50, 50] %[xmin ymin xmax ymax]


recon = ["AI","IR", "FBP"];
dia = [16,21,26,31,36];
psize = [0.7832, 0.7832];             % pixel 
ny = 1./(2*psize);              % Nyquist frequency

% n = length(Files);
% for i=1:n
%     File = strcat(Path,Files(i));
%     I = dicomread(File); 
%     [f,nps] = npsplot(I, ROI,psize);
%     Q = cumtrapz(f,nps);
%     X(:,i) = f;
%     Y(:,i) = smooth(Q,10);
% %     leg{i} = [num2str(AI_set(i)) '%'];
%     leg{i} = [strcat(recon(i))];
% end
k = length(dia);
n = length(recon);
for j=1:k 
    for i=1:n
        File = strcat(Path,Files30(j,i));
        I = dicomread(File);
        [f,nps] = npsplot(I,ROI,psize);
        Q_30(j,i) = trapz(nps);
        T = f .* nps;
        CFR_30(j,i) = trapz(T)/Q_30(j,i);
        X(j,i) =dia(j);
    end

end

for j=1:k 
    for i=1:n
        File = strcat(Path,Files60(j,i));
        I = dicomread(File);
        [f,nps] = npsplot(I,ROI,psize);
        Q_60(j,i) = trapz(nps);
        T = f .* nps;
        CFR_60(j,i) = trapz(T)/Q_60(j,i);
        X(j,i) =dia(j);
    end

end

for j=1:k 
    for i=1:n
        File = strcat(Path,Files90(j,i));
        I = dicomread(File);
        [f,nps] = npsplot(I,ROI,psize);
        Q_90(j,i) = trapz(nps);
        T = f .* nps;
        CFR_90(j,i) = trapz(T)/Q_90(j,i);
        X(j,i) =dia(j);
    end

end
CFR_30 = CFR_30(:,any(CFR_30));
CFR_60 = CFR_60(:,any(CFR_60));
CFR_90 = CFR_90(:,any(CFR_90));
% to test validity

k = length(dia);
n = length(recon);
for j=1:k 
    
   A(:,j) = sqrt(Q_30(j,1)/Q_30(j,2))
   B(:,j) = Q_60(j,1)/Q_60(j,2)
   C(:,j) = Q_90(j,1)/Q_90(j,2)

end



% loop for CFR
k = 2;
for i=1:k
    CFR_3(:,i) = CFR_30(:,i)./CFR_30(:,3);
    CFR_6(:,i) = CFR_60(:,i)./CFR_60(:,3);
    CFR_9(:,i) = CFR_90(:,i)./CFR_90(:,3);
end
CFR_3 = CFR_3(:,any(CFR_3));
%CFR_3(:,3) = [];
CFR_6 = CFR_6(:,any(CFR_6));
%CFR_6(:,3) = [];
CFR_9 = CFR_9(:,any(CFR_9));
%CFR_9(:,3) = [];
% loop for NMR
k = 2;
for i=1:k
    NMR_30(:,i) = sqrt(Q_30(:,i) ./ Q_30(:,3));
    NMR_60(:,i) = sqrt(Q_60(:,i) ./ Q_60(:,3));
    NMR_90(:,i) = sqrt(Q_90(:,i) ./ Q_90(:,3));  
end
NMR_30 = NMR_30(:,any(NMR_30));
NMR_60 = NMR_60(:,any(NMR_60));
NMR_90 = NMR_90(:,any(NMR_90));

 % show ROI on slice with smallest diameter
%  K_file = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020018'
%  K = dicomread(K_file);
%  imshow(K,[])
%  hold on;
%  rectangle('Position',ROI_display, 'EdgeColor', 'r', 'LineWidth', 2);
%  hold off;

%plot(X,Y,'LineWidth',1.5);
% trying to plot NMR on y axis and CFR on X axis like the paper
% will have 30%, 60%, and 90% AI and IR 
CFR = [CFR_3, CFR_6, CFR_9];
NMR = [NMR_30, NMR_60, NMR_90];
labels = ["30% AI" "30% IR" "60% AI" "60% IR", "90% AI" "90% IR"];
% scatter(dia,CFR);
% scatter(NMR_30, CFR_3)
% scatter(dia, NMR);
% scatter(dia,CFR, 'filled', 'LineWidth',1.5);
markerSize = 100;
% scatter(dia, CFR, markerSize, 'filled');
% scatter(dia, NMR, markerSize, 'filled');
scatter(CFR,NMR, markerSize, 'filled');
legend(labels);
ylabel('NMR');
xlabel('CFR');
title('NMR and CFR');

function [f, nps, noise, fpeak, fav] = npsplot(I,ROIs,psize,varargin)
%This function measures the 2D NPS from an image

%ROIs [xmin ymin xsize ysize]

%get the inputs
[padSize,subtractMethod,fRange,fSamples] = parseInputs(varargin);

%get the number of slices
slices = size(I,3);

%create the frequency axes
fx = getFFTfrequency(psize(1),padSize);
fy = getFFTfrequency(psize(2),padSize);
[Fx,Fy]=meshgrid(fx,fy);
[~,Fr]=cart2pol(Fx,Fy);
dfx = fx(2)-fx(1);
dfy = fy(2)-fy(1);

%loop over all ROIs and slices
n=0;
i = 1;
ROI = ROIs(i,:);
N = prod(psize)/prod(ROI(3:4)); %normalization factor for this ROI
for j=1:slices
    n=n+1;
    %get the pixels
    im = I(ROI(2):ROI(4),ROI(1):ROI(3),j);
    %imshow(im)
    %remove the mean
    im = subtractMean(im,subtractMethod,psize);
        
    NPS(:,:,n) = N*abs(fftshift(fft2(im,padSize,padSize))).^2;
        
    nps=NPS(:,:,n);
    variance(n,1) = sum(nps(:)*dfx*dfy);
        
end


%Average over all ROI positions
NPS = mean(NPS,3);

%compute the variance
noiseSTD = std(sqrt(variance));
variance=mean(variance);
noise = sqrt(variance);

%Get radial average
edges = linspace(fRange(1),fRange(2),fSamples);
[f, nps, rvar] = rebinData(Fr, NPS, edges, 0);
f=f'; nps=nps';

%get the peak frequency
[~,ind] = max(smooth(nps));
fpeak = f(ind);     %Peak of 1D radial NPS

%Get the average frequency from 1D radial NPS
p = nps/sum(nps);   %turn 1D nps into probability distribution
fav = sum(p.*f); 

end



function [padSize,subtractMethod,fRange,fSamples] = parseInputs(args)

padSize = 512;
subtractMethod = 'poly'; %'mean';
fRange = [0 2];
fSamples = 128;
if ~isempty(args)
    for i=1:2:length(args)
        switch args{i}
            case 'padSize'
                padSize = args{i+1};
            case 'subtractMethod'
                subtractMethod = args{i+1};
            case 'fRange'
                fRange = args{i+1};
            case 'fSamples'
                fSamples = args{i+1};
        end
    end
end
end
