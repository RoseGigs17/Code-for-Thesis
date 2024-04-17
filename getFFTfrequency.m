% Copyright 20xx - 2019. Duke University
function f = getFFTfrequency(psize,N,varargin)

switch length(varargin)
    case 1
        mode = varargin{1};
    otherwise
        mode = 'shifted'; %Default case
        
end
    Fs = 1/psize;   %Sampling rate
    df = Fs/N;      %spacing in frequency domain
    n=0:(N-1);      %Vector of indices
    f = n*df;       %frequency (unshifted version)
    
    switch mode
        case 'shifted'
            fnew = fftshift(f); %get the shifted version of f
            ind = fnew==0; %find index of dc component in shifted FFT
            f = f - f(ind); %Shift f such that 0 shows up at dc bin component 
    end

end