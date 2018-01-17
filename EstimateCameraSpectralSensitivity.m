close all;
% clear all;

%% Recover camera spectral sensitivity from 20 LED images
% What is needed:
% - Principal component database
% - Dark image
% - 20 images of single LED channels
% - 20 spectra of single LED channels
% - White image (D65)
% - Spectrum of white illuminant (D65)
% - (Region of interest (ROI) for each image)

%% Principal component database

% Get the eigenvectors of the camera spectral sensitivity
numEV = 3;

% [eRed,eGreen,eBlue]=PCACameraSensitivity(numEV);

freq_range = [380; 821];
[freqPC, eRed, eGreen, eBlue] = GeneratePCAModel('testdata\responses', numEV, freq_range);

% freq_range = [400; 720];
% [freqPC, eRed, eGreen, eBlue] = loadPCADB('testdata\cameraPCs_correctOrder.txt');


%% Load images and spectra

% Dark image
% The dark image is important for compensation of potential offset in digital value range and (fixed pattern)
% noise. Although the analyzer class does not enforce to use a dark image, it is strongly advised.
% The dark image requires a ROI, you may use the following code for placing a ROI in the center or create
% your own ROI object:
darkImage = imread('testdata\spectral-led_0_0_25000.pgm');
[h, w, c] = size(darkImage);

% White image (D65) + Spectrum of white illuminant (D65)
% White images are needed for a correct white balancing of the spectral sensitivities of a color camera.
% Since the spectral sensitivities for each color channel are estimated separately, at least one white
% spectrum is needed for the determination of the channel ratios in the algorithm of the spectral sensitivity
% analysis.
whiteImage = imread('testdata\spectral-D65_200_25000.pgm');
whiteImage = whiteImage - darkImage;
whiteImageRGB = zeros(h/2, w/2, 3);
g1 = whiteImage(1:2:end, 1:2:end);
b = whiteImage(1:2:end, 2:2:end);
r = whiteImage(2:2:end, 1:2:end);
g2 = whiteImage(2:2:end, 2:2:end);
whiteImageRGB(:,:,1) = r;
whiteImageRGB(:,:,2) = (g1 + g2) / 2;
whiteImageRGB(:,:,3) = b;

[freq, whiteSpectrum] = loadSpectrum('testdata\D65.measured');
N = 1 + freq_range(2) - freq_range(1);
whiteSpectrum = InterpolateSpectraTo1nm(whiteSpectrum, [freq(1); freq(end); freq(2)-freq(1)], freq_range);

% 20 images of single LED channels + spectra
imageLedRGB = zeros(20, h/2, w/2, 3);
spectrumLED = zeros(20, N);
for i=1:20
    imageFilename = strcat('testdata\spectral-led_', num2str(i), '_1000_25000.pgm');
    spectrumFilename = strcat('testdata\led', num2str(i), '_1000.000000.measured');
    img = imread(imageFilename);
    [f, spectrum] = loadSpectrum(spectrumFilename);
    
    if isequal(freq, f) ~= 1
       error(strcat('Spectrum for file ', spectrumFilename, ' has wrong range!'));
    end
    
    % Subtract dark image here? Yes! 
    img = img - darkImage;
    g1 = img(1:2:end, 1:2:end);
    b = img(1:2:end, 2:2:end);
    r = img(2:2:end, 1:2:end);
    g2 = img(2:2:end, 2:2:end);
    imageLedRGB(i, :, :, 1) = r;
    imageLedRGB(i, :, :, 2) = (g1 + g2) / 2;
    imageLedRGB(i, :, :, 3) = b;
    spectrumLED(i, :) = InterpolateSpectraTo1nm(spectrum, [freq(1); freq(end); freq(2)-freq(1)], freq_range);
end

% Normalize image values to [0-1] range
imageLedRGB = imageLedRGB / (2^16 - 1);

figure,
w = freq_range(1):freq_range(2);
hold on
for i=1:20
    plot(w,spectrumLED(i,:));
    hold on;
end
title('LED spectra');
hold off

%% Set ROI and get mean value for each LED image
roiW = 100;
roiH = 100;
imageLedRgbRoiMean = zeros(20, 3);
startY = h/4 - roiH/2;
startX = w/4 - roiW/2;
for i=1:20
    roi = imageLedRGB(i,startY:startY+roiH-1, startX:startX+roiW-1, :);
    r = squeeze(roi(:,:,:,1));
    g = squeeze(roi(:,:,:,2));
    b = squeeze(roi(:,:,:,3));
    meanR = mean(r(:));
    meanG = mean(g(:));
    meanB = mean(b(:));
    imageLedRgbRoiMean(i, :) = [meanR meanG meanB];
end

% Same for white image
roi = whiteImageRGB(startY:startY+roiH-1, startX:startX+roiW-1, :);
r = roi(:,:,1);
g = roi(:,:,2);
b = roi(:,:,3);
meanR = mean(r(:));
meanG = mean(g(:));
meanB = mean(b(:));
whiteImageRgbRoiMean = [meanR meanG meanB];

%% Recover camera spectral sensitivity from a single image under known daylight
clear cmfHat

[cmfHat(:,1)] = RecoverCSSwithLED(whiteImageRGB(:,:,1), whiteSpectrum, imageLedRgbRoiMean(:,1), spectrumLED, eRed);
[cmfHat(:,2)] = RecoverCSSwithLED(whiteImageRGB(:,:,2), whiteSpectrum, imageLedRgbRoiMean(:,2), spectrumLED, eGreen);
[cmfHat(:,3)] = RecoverCSSwithLED(whiteImageRGB(:,:,3), whiteSpectrum, imageLedRgbRoiMean(:,3), spectrumLED, eBlue);

cmfHat=cmfHat./max(cmfHat(:));
cmfHat(cmfHat<0)=0;

% Calculate channel ratios using D65 spectrum
rXYZ = sum(cmfHat(:,1) .* whiteSpectrum);
gXYZ = sum(cmfHat(:,2) .* whiteSpectrum);
bXYZ = sum(cmfHat(:,3) .* whiteSpectrum);

rGain = meanR / rXYZ;
gGain = meanG / gXYZ;
bGain = meanB / bXYZ;

cmfHatWB(:,1) = cmfHatWB(:,1) * rGain;
cmfHatWB(:,2) = cmfHatWB(:,2) * gGain;
cmfHatWB(:,3) = cmfHatWB(:,3) * bGain;

cmfHatWB=cmfHatWB./max(cmfHatWB(:));
cmfHatWB(cmfHatWB<0)=0;

%% Load the measured cmf (ground truth) of the camera
% camName='Canon60D';
% [rgbCMF,camNameAll]=getCameraSpectralSensitivity();
% for i=1:length(camNameAll)
%     if(strcmp(camNameAll{i},camName))
%         cmf=[rgbCMF{1}(:,i),rgbCMF{2}(:,i),rgbCMF{3}(:,i)];
%     end
%     
% end
% cmf=cmf./max(cmf(:));
% w=400:10:720;

% Read CSV, assuming GRBG order
clear cmf;
data = csvread('testdata\responses\center_cam0_correctExposure.csv');
f = data(:, 1);
cssR = data(:, 3);
cmf(:,1) = InterpolateSpectraTo1nm(cssR, [f(1); f(end); f(2)-f(1)], freq_range);
cssG = (data(:, 2) + data(:, 5)) / 2;
cmf(:,2) = InterpolateSpectraTo1nm(cssG, [f(1); f(end); f(2)-f(1)], freq_range);
cssB = data(:, 4);
cmf(:,3) = InterpolateSpectraTo1nm(cssB, [f(1); f(end); f(2)-f(1)], freq_range);
w = freq_range(1):freq_range(2);

cmf=cmf./max(cmf(:));

%% Plot estimated response vs ground truth

figure;
plot(w,cmf(:,1),'r');
hold on;
plot(w,cmf(:,2),'g');
hold on;
plot(w,cmf(:,3),'b');
hold on;

plot(w,cmfHat(:,1),'r--');
hold on;
plot(w,cmfHat(:,2),'g--');
hold on;
plot(w,cmfHat(:,3),'b--');
hold off;

title('Original camera response functions vs estimated');
legend('R','G','B','R_e','G_e','B_e');


figure;
plot(w,cmf(:,1),'r');
hold on;
plot(w,cmf(:,2),'g');
hold on;
plot(w,cmf(:,3),'b');
hold on;

plot(w,cmfHatWB(:,1),'r--');
hold on;
plot(w,cmfHatWB(:,2),'g--');
hold on;
plot(w,cmfHatWB(:,3),'b--');
hold off;

title('Original camera response functions vs estimated (WB)');
legend('R','G','B','R','G_b','B_b');

% legend('r_m','g_m','b_m','r_e','g_e','b_e');
