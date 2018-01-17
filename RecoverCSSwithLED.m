
%Summary
%   The function is to recover the camera spectral sensitivity given (the
%   spectral reflectance of the samples,) the eigenvectors of the camera
%   sensitivity, and the illuminant spectrum
%
%[IN]
%   ill: the light source spectrum
%   reflSet: the spectral reflectance of samples
%   w: wavelength range
%   XYZSet: the radiance captured by the camera
%   e: eigenvector of the camera spectral sensitivity
%
%[OUT]
%   X: the recovered camera spectral sensitivity
%   A, b: Ax=b
%
function [X,A,b]=RecoverCSSwithLED(whiteImageCh, whiteSpectrum, roisCh, spectrumLED, eCh)


numLEDs = size(roisCh, 1);
A = zeros(numLEDs, size(eCh,2));  % size(eCh,2) is 3 because of 3 principal components (and not because of 3 channels!)
b = zeros(size(A,1), 1);

% deltaLambda=10;

for i=1:numLEDs
    a1 = spectrumLED(i,:); %diag(spectrumLED(i,:));
    a2 = eCh;
    b1 = roisCh(i);
    
    A(i,:) = spectrumLED(i,:) * eCh;% .* deltaLambda;
    b(i) = roisCh(i);
end

% numRefl=size(reflSet,2);
% 
% A=zeros(numRefl,size(e,2));
% b=zeros(size(A,1),1);
% 
% deltaLambda=10;
% 
% for i=1:numRefl
%     A(i,:)=reflSet(:,i)' * diag(ill) * e .* deltaLambda;
%     b(i)=XYZSet(i);
% end


% I = C * L * R; % pixel_value = camera_spectral_sensitivity * illuminant * reflectance
% I = sigma * E * L * R;
% A = ELR; % eigenvectors * illuminant * reflectance
% b = ik;

X = A \ b;

X = eCh * X;


end
