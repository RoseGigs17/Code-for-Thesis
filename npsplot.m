% Example Input
   File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010002';
   I = dicomread(File); 
   ROI = [100,100,400,400]
   psize = [0.1, 0.1]
%
% Call:  [f,nps] = npsplot(I, ROI,psize);

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

padSize = 128;
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
