function result = InterpolateSpectraTo1nm(input_values, input_range, output_range)
%Interpolates input spectra to 1nm interval
%
%input_values = Spectra, Spectrum in each column
%
%input_range = [1st input sample waveleght;
%               last input sample waveleght;
%               sampling_interval]
%or same data in columns for all spectra in columns
%
%output_range = [1st output sample waveleght;
%                last output sample waveleght]
%
%result = 1nm sampled spectra, Spectrum in each column
%
%If input_range has only one column same sampling and range is
%used for all input spectra
%
%Output range needs to be defined as integer values (1nm accuracy)

%If only one column in input range, repeat for all columns
if(size(input_range,2) == 1)
    input_range_colums = repmat(input_range, 1, size(input_values,2));
    
%else use input_range directly
%(no validity check, must have correct size)
else
    input_range_colums = input_range;
end

%Initialize result with zeros
result = zeros((1 + output_range(2,:)-output_range(1,:)),size(input_values,2));

%Result is sampled with 1nm interval, 1st value is the start wavelength,
%last value is the end wavelenght
result_sampling = output_range(1,:):output_range(2,:);

%Loop through all spectra in columns
for spectrum_index = 1:size(input_values,2)
    
    %If there is some input range given
    if(input_range_colums(2,spectrum_index)-input_range_colums(1,spectrum_index) > 0),
        
        %Input sampling
        input_sampling = input_range_colums(1,spectrum_index): ...
            input_range_colums(3,spectrum_index): ...
            input_range_colums(2,spectrum_index);
        
        %Interpolate result to 1nm interval, fill values outside input range with 0
        %Method can be chosen, initially 'linear'
        result(:,spectrum_index) = ...
            interp1(input_sampling, ...
                    input_values(1:length(input_sampling), ...
                                 spectrum_index), ...
                    result_sampling, 'linear', 0);
        
        %If input range was 0 or less
    else        
        %result column is full of zeros
        result(:,spectrum_index) = zeros((1 + output_range(2,:)-output_range(1,:)),1);
    end
end

end