% 0-percent or FBP PATH
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010002/';
% slices =["D00020018" "D00020043" "D00020069" "D00020080" "D00020109"];

% 30-percent AI PATH
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010003/';
% slices =["D00030018" "D00030043" "D00030069" "D00030080" "D00030109"];

% 60-percent AI PATH
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010004/';
% slices =["D00040018" "D00040043" "D00040069" "D00040080" "D00040109"];

% 90-percent AI PATH
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010005/';
% slices =["D00050018" "D00050043" "D00050069" "D00050080" "D00050109"];

% % 30-percent IR PATH
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/';
% slices =["D0018" "D0043" "D0069" "D0080" "D0109"];

% 60-percent IR PATH
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/';
% slices =["D0018" "D0043" "D0069" "D0080" "D0109"];

% 90-percent IR PATH
Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/';
slices =["D0018" "D0043" "D0069" "D0080" "D0109"];


%large ROI
% ROI = [200,210,320,330]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
% ROI_display = [200, 210, 120, 120] %[xmin ymin xmax ymax]

%small ROI no middle rod
ROI = [270,266,320,316]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
ROI_display = [270, 266, 50, 50] %[xmin ymin xmax ymax]

% large ROI no middle rod
% ROI = [190,280,330,330]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
% ROI_display = [190, 280, 140, 50] %[xmin ymin xmax ymax]

%bone ROI
% ROI = [310,240,340,250];
% ROI_disp = [310,240,30,10];

%small ROI middle rod
% ROI = [250,266,300,316]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
% ROI_display = [250, 266, 50, 50] %[xmin ymin xmax ymax]

Dia = [16, 21, 26, 31, 36];
psize = [0.7832, 0.7832];             % pixel 
ny = 1./(2*psize);              % Nyquist frequency

n = length(slices);
for i=1:n
    File = strcat(Path,slices(i));
    I = dicomread(File); 
    [f,nps, noise] = npsplot(I, ROI,psize);
    X(:,i) = f;
    Y(:,i) = smooth(nps, 10);
    %Y(:,i) = nps;
    leg{i} = [num2str(Dia(i)) ' cm'];
end

 % show ROI on slice with smallest diameter
%  K_file = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020018'
%  K = dicomread(K_file);
%  imshow(K,[-500,-50])
%  hold on;
%  rectangle('Position',ROI_display, 'EdgeColor', 'r', 'LineWidth', 2);
%  hold off;

%visual the difference in noise
% S_file = strcat(Path,slices(5));
% S = dicomread(S_file);
% S_detail = S(ROI(2):ROI(4),ROI(1):ROI(3));
% imshow(S)
%imshow(S_detail)

plot(X,Y,'LineWidth',1.5);
legend(leg);
ylabel('NPS [mm^2HU^2]');
xlabel('Spatial Frequency [mm^-^1]');
title('NPS - 90% IR');
xlim([0 1])

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
