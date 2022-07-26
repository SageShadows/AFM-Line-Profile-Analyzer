function [thicknesses, areas] = Area_AFMProfileAnalyzer(varargin)
% 
% Inputs:
% 
% isFile: '' if you're inputting a file
%         [A] if you want to use data in [thicknesses, lengths] format that
%         is formatted properly (both line profiles in sequential manner.
% 
% Outputs:
% thicknesses: a vector containing thicknesses of accepted flakes in nm
%              averaged over both values provided
% 
% areas: a vector containing areas of accepted flakes in nm^2. 

% v1.0: First implementation of code as simple data importer and script 
% analyzer. 2/12/2019
% David Lam, (c) 2019

%% Data Import 
% Prompts user to select files of interest.
    data = [];
    if isempty(varargin{1})
        [file, path] = uigetfile('', 'Select a file'); 
        fileID = fopen(strcat(path, file));
        currentLine = fgetl(fileID);
        currentLine = fgetl(fileID);
        data = [];
        while currentLine ~= -1
            currentRow = textscan(currentLine, '%s');
            data = [data; str2double(currentRow{1,1})'];
            currentLine = fgetl(fileID);
        end
        fclose(fileID);
    else
        data = varargin{1}; 
    end
    
%% Data Analysis
% Analyzes data by first reshaping the original data matrix, then doing
% consistency checks + multiplies the axes lengths for a rough estimate of
% the area. 
    data = reshape(data', 4, size(data,1)/2)';
    thicknesses = [];
    areas = [];
    for i = 1:size(data,1)
        userAnswer = questdlg(sprintf('Thickness 1: %3.2f nm, Thickness 2: %3.2f nm, Delta T: %3.2f nm'...
            , [data(i,1), data(i,3), abs(data(i,1)-data(i,3))]));
        drawnow; %Stops UI lag from the questdlg
        close
        if userAnswer == "Yes"
            thicknesses = [thicknesses mean([data(i,1), data(i,3)])];
            areas = [areas areaFunction(data(i,2), data(i,4))];
        elseif userAnswer == "No"
            continue
        elseif userAnswer == "Cancel"
            break
        end
    end
    
    function [flakeArea] = areaFunction(a1, a2, model)
        flakeArea = pi*a1*a2/4; %Elliptical model 
    end
end

% Outputs a histogram for quick viewing. Uncomment if you want this feature.  
%{
close all; 
figure(1); histogram(thicknesses); title("Thickness Histogram"); 
xlabel('Thickness [nm]'); ylabel('Frequency');

figure(2); histogram(lengths); title("Length Histogram");
xlabel('Lateral Length [nm]'); ylabel('Frequency');
%}
