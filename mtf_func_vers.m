% --- Load in Files --- %

%load in 16 cm - 0% or FBP 
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020022';
% I = dicomread(File);

%load in 16 cm - 30% AI 
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010003/D00030022';
% I = dicomread(File);

%load in 16 cm - 60% AI 
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010004/D00040022';
% I = dicomread(File);

%load in 16 cm - 90% AI 
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/00010005/D00050022';
% I = dicomread(File);

%load in 16 cm - 30% IR 
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0022';
% I = dicomread(File);

%load in 16 cm - 60% IR 
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0022';
% I = dicomread(File);

%load in 16 cm - 90% IR 
% File = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0022';
%I = dicomread(File);

%try creating path to asses all levels of intensity for one method of
%reconstruction

%AI
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/Mercury4/Neusoft AI Mercury4 Image Data/';
% levels = ["00010002/D00020022","00010003/D00030022","00010004/D00040022","00010005/D00050022"];
%IR
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
% levels = ["Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020022","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0022","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0022","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0022"];

% to compare alike levels AI first, IR Second

% 30%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
% levels = ["Mercury4/Neusoft AI Mercury4 Image Data/00010003/D00030022","ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0006/D0022", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020022"];

% 60%
% Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
% levels = ["Mercury4/Neusoft AI Mercury4 Image Data/00010004/D00040022", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0007/D0022", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020022"];

% 90%
Path = '/Users/rosegigliotti/Documents/Grad School Stuff/Research/RESEARCH/';
levels = ["Mercury4/Neusoft AI Mercury4 Image Data/00010005/D00050022", "ClearView_Iterative_Recon_Test/MercuryPhantomclearview/PAT000_2/S0008/D0022", "Mercury4/Neusoft AI Mercury4 Image Data/00010002/D00020022"];

% --- ROI Selection --- %

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
ROI = [250,207,280,217];
ROI_disp = [250,207,30,10];

% -- legend stuff -- %
%if looking at all intensity
% legend_key = [0,30,60,90];

%if comparing ai to ir of same level
legend_key = ["AI" "IR" "FBP or baseline"]

n = length(levels);

for i=1:n
    File = strcat(Path,levels(i))
    [freq,mtf_y] = mtf(File,ROI);
    X(:,i) = freq;
    Y(:,i) = smooth(mtf_y,10);
    %if numerical use this
    leg{i} = [num2str(legend_key(i)) '%'];
    %if string use this
%     leg{i} = [strcat(legend_key(i))];
end

%graph

plot(X,Y,'LineWidth',1.5);
legend(leg);
ylabel('\fontsize{12}MTF [%]');
xlabel('\fontsize{12}Spatial Frequency [lp/mm]');
title('\fontsize{16}MTF - 90% - Solid Water Insert');
%xlim([0 1])



function[freq, mtf_y] = mtf(File,ROI)
%gives you the mtf by taking an roi and going from esf to lsf to mtf

%ESF part
I = dicomread(File);
B = I(ROI(2):ROI(4),ROI(1):ROI(3));
px_size = 0.7832031; %mm
mean_B = mean(B);
B_col = mean_B';
x_col = (1:length(B_col))';
x_dist = x_col * px_size; %convert px size to mm

%LST part
lst_y_raw = (diff(B_col));
lst_y = abs((diff(B_col)));
lst_x = (1:30)' * px_size;
%plot(lst_x,lst_y);
lst = [lst_x lst_y];

%MTF Part
FFTf_x=fft((lst_y));
MTF=abs(real(FFTf_x(1:end))/real(FFTf_x(1)));
freq = [0:(length(lst_y)-1)/2]/10.;
mtf_y = MTF(1:length(freq));
end


