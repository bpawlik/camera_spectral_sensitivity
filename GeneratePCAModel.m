
%Summary
%   The function is to do PCA on camera sensitivity
%[IN]
%   numEV: number of eigenvectors to retain
%
%[OUT]
%   eRed, eGreen, eBlue: the eigenvectors of each of the three channel
%
function [freq, eRed,eGreen,eBlue]=GeneratePCAModel(folder, numEV, freq_range)

% Read all camera spectral responses from the folder
myFiles = dir(fullfile(folder,'*.csv'));
N = length(myFiles);
freq = freq_range(1):freq_range(2);
freqN = length(freq);

redCMF = zeros(freqN, N);
greenCMF = zeros(freqN, N);
blueCMF = zeros(freqN, N);

for k = 1:N
    baseFileName = myFiles(k).name;
    fullFileName = fullfile(folder, baseFileName);
    
    % Read CSV, assuming GRBG order
    data = csvread(fullFileName);
    f = data(:, 1);
    cssR = data(:, 3);
    cssRinterp1nm = InterpolateSpectraTo1nm(cssR, [f(1); f(end); f(2)-f(1)], freq_range);
    cssG = (data(:, 2) + data(:, 5)) / 2;
    cssGinterp1nm = InterpolateSpectraTo1nm(cssG, [f(1); f(end); f(2)-f(1)], freq_range);
    cssB = data(:, 4);
    cssBinterp1nm = InterpolateSpectraTo1nm(cssB, [f(1); f(end); f(2)-f(1)], freq_range);
    
    redCMF(:, k) = cssRinterp1nm;
    greenCMF(:, k) = cssGinterp1nm;
    blueCMF(:, k) = cssBinterp1nm;
end

figure,
hold on
for i=1:size(greenCMF,2)
    plot(freq,redCMF(:,i),'r');
    hold on;
    plot(freq,greenCMF(:,i),'g');
    hold on;
    plot(freq,blueCMF(:,i),'b');
    hold on;
end
title('Original camera response functions (whole set)');
hold off

%normalize to each curve
for i=1:size(greenCMF,2)
    redCMF(:,i)=redCMF(:,i)./max(redCMF(:,i));
    greenCMF(:,i)=greenCMF(:,i)./max(greenCMF(:,i));
    blueCMF(:,i)=blueCMF(:,i)./max(blueCMF(:,i));
end


%% do PCA on cmf

if(nargin>0)
    retainEV=numEV;
else
    retainEV=1;
end


[eRed]=GetEigenvector(redCMF,retainEV);
[eGreen]=GetEigenvector(greenCMF,retainEV);
[eBlue]=GetEigenvector(blueCMF,retainEV);


end
