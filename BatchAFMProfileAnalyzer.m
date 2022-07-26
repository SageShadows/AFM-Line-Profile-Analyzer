function [vertDim, latDim] = BatchAFMProfileAnalyzer(programMode, varargin)
% BatchAFMProfileAnalyzer is a wrapper function that simplifies operation 
% of AFMProfileAnalyzer by introducing a UI interface to select multiple
% files for analysis. Read AFMProfileAnalyzer for more in-depth
% documentation. 
% Inputs:
% programMode: a string that indicates operation mode. Use one of the below.
% "auto": runs code without supervision. Not advised. 
% "noErrors": does not add data points with errors outputted in the command
%             window. OK for rapid analysis. 
% "manual" : views each profile and allows you to accept or reject for
%           inclusion. In this mode, a plot of the profile and the second 
%           derivative are presented. The dashed lines indicate where the 
%           cutoff for the flake will be. Best mode to use for quality
%           control. 
% "area" : runs code in manual mode, but pairs off the profiles as they are
% given in major/minor axis form. 
%
% outputFile: an optional string that indicates desired output file. If no
%             output file is specified, will simply output the data into
%             MATLAB instead of saving to a file. 
% 
% Outputs:
% thicknesses: a vector containing thicknesses of accepted flakes in nm. 
% 
% lengths: a vector containing lengths of accepted flakes in nm. 

% v1.2: Updated to include "area" as programMode. Functionality is
% abstracted to AFMProfileAnalyzer.m; edits to this wrapper function are
% cosmetic in nature. 
% v1.1: Minor fixes to code, inclusion of histogram feature to output
% histograms (default commented out) 5/23/2018
% v1.0: First implementation of code. 5/16/2018
% David Lam, (c) 2018

% Initializes output information. 
vertDim = [];
latDim = [];

% Prompts user to select files of interest.
[file, path] = uigetfile('.txt', 'Select One or More Files', 'MultiSelect', 'on'); 

% If more than one file selected, iterates through each file. 
if iscell(file)
    for i = 1:length(file)
        [v, l] = AFMProfileAnalyzer(strcat(path, file{i}),programMode);
        vertDim = [vertDim; v];
        latDim = [latDim; l]; 
    end
else 
    [vertDim, latDim] = AFMProfileAnalyzer(strcat(path, file),programMode);
end

% Writes output to desired file if file name specified. 
if nargin == 2
    s1 = "Thickness [nm]";
    if programMode == "area"
        s2 = "Area [nm^2]";
    else
        s2 = "Lateral Length [nm]";
    end
    fileID = fopen(varargin{1},'w');
    fprintf(fileID,'%s \t %s \n', [s1; s2]);
    fprintf(fileID,'%3.2f \t %4.2f\n', [vertDim'; latDim']);
    fclose(fileID);
end 

% Outputs a histogram for quick viewing. Uncomment if you want this feature.  
%{
close all; 
figure(1); histogram(thicknesses); title("Thickness Histogram"); 
xlabel('Thickness [nm]'); ylabel('Frequency');

figure(2); histogram(lengths); title("Length Histogram");
xlabel('Lateral Length [nm]'); ylabel('Frequency');
%}
