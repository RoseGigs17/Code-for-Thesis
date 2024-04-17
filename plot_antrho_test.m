% Anthro File - 0% or FBP
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/AI Testing/S0002/D0056';
% Anthro File - 30% AI
%File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/AI Testing/S0003/D0056';
% Anthro File - 60% AI
%File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/AI Testing/S0004/D0056';
% Anthro File - 90% AI
%File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/AI Testing/S0005/D0056';

Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/AI Testing/';
Files = ["S0002/D0056" "S0003/D0056" "S0004/D0056" "S0005/D0056"];


%ROI for Anthro - Homogenous
ROI = [120,120,180,180]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
ROI_display = [120, 120, 60, 60] %[xmin ymin xmax ymax]


%ROI for antrho - heterogenous 
% ROI = [295,285,345,335]; %[xmin ymin (xmin+xmax) (ymin+ymax)]
% ROI_display = [295, 285, 50, 50] %[xmin ymin xmax ymax]

AI_set = [0, 30, 60, 90];
dia = [16,21,26,31,36]
psize = [0.7832, 0.7832];             % pixel 
ny = 1./(2*psize);              % Nyquist frequency

n = length(Files);
for i=1:n
    File = strcat(Path,Files(i));
    I = dicomread(File); 
    [f,nps] = npsplot(I, ROI,psize);
    X(:,i) = f;
    Y(:,i) = smooth(nps,10);
    leg{i} = [num2str(AI_set(i)) '%'];
end

 % show ROI on slice with smallest diameter
%  K_file = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/AI Testing/S0002/D0056'
%  K = dicomread(K_file);
%  imshow(K,[])
%  hold on;
%  rectangle('Position',ROI_display, 'EdgeColor', 'r', 'LineWidth', 2);
%  hold off;

plot(X,Y,'LineWidth',1.5);
legend(leg);
ylabel('NPS [mm^2HU^2]');
xlabel('Spatial Frequency [mm^-^1]');
title('Noise Power Spectra - Homogenous Region');
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
