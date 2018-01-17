function [frequency, spectrum] = loadSpectrum(filename)

data = dlmread(filename, ' ', [8 0 449 1]);
frequency = data(:,1); % column 1 of the data text file is assigned the variable x
spectrum = data(:,2); % column 2 is assigned the variable y

end