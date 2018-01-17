function [frequency, eRed, eGreen, eBlue] = loadPCADB(filename)

data = dlmread(filename, ' ', [11 0 331 9]);
frequency = data(:,1); % column 1 of the data text file is assigned the variable x
eRed = data(:,2:4); % column 2 is assigned the variable y
eGreen = data(:,5:7);
eBlue = data(:,8:10);
end