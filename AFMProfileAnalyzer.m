function [vertDim, latDim] = AFMProfileAnalyzer(fileName,programMode)
% AFMProfileAnalyzer is meant to analyze large numbers of AFM flake
% profiles.
% 
% Inputs:
% fileName: a string indicating the input text file. The formatting for an 
%           input text file follows the order of multiple profiles  
%           concatenated in the following structure: 
%           [xdir_1 ydir_1 ... xdir_final ydir_final]  
%           and is robust enough to handle profiles of variable points  
%           so long as the length of xdir_j and ydir_j are the same.
% 
% programMode: a string that indcates operation mode. Use one of the below.
% "auto": runs code without supervision. Not advised. 
% "noErrors": does not add data points with errors outputted in the command
%             window
% "manual" : views each profile and allows you to accept or reject for
%           inclusion. In this mode, a plot of the profile and the second 
%           derivative are presented. The dashed lines indicate where the 
%           cutoff for the flake will be. Best mode to use for quality
%           control. 
% "area" : runs code in manual mode, but pairs off the profiles as they are
% given in major/minor axis form. 
% 
% Outputs:
% vertDim: a vector containing thicknesses of accepted flakes in nm. 
% 
% latDim: a vector containing lengths of accepted flakes in nm or areas of
% accepted flakes in nm^2.

% Ideally, input files are generated in Gwyddion by drawing line profiles
% over single flakes of interest via the "Extract Profiles" option
% and generating a single plot of all of the profiles by unselecting the
% "Separate profiles" option. Then, click apply to generate the profile plot
% and right click on the profile. Select "Export Text", unselect all of the
% options such that you only have raw data, then click OK and save the file.

% At this time, only individual profiles of a single flake can be
% analyzed. To find the area of a flake, you must input a file with pairs
% of profiles (i.e. profiles 1 and 2 are of the same flake in orthogonal
% directions)

% DO NOT IMPORT PROFILES WITH TWO OR MORE FLAKES!!
% The analysis works as follows: It determines the global maximum in the
% profile, then scans in a range around the global maximum to find the local
% minimum for the second derivative. Then, the code determines the second
% derivative maximum to the left and right of this calculated local second
% derivative minimum. Afterwards, the flake cutoff is determined by finding
% the position that is a certain percentage of the second derivative maxima
% on the left and right within some range of the maxima. 

% Operation of this script is as follows: 
% Choose a file that follows the same structure as above. The script will
% extract the text and form a matrix of the profile information. Then,
% choose the program mode. Options available are outlined in the "User
% Imports" section following this intro section. Other parameters should not 
% be changed unless you understand the code, but those parameters can be 
% found under the "Profile Analysis" section of this script. Recommended
% operation mode for this program is in "manual" mode, which allows for
% some measure of user input. 

% v1.2: Included an area function that models a flake as an ellipse, where
% the profiles inputting are the major and minor axes of an ellipse.
% 2/21/2019
% v1.1: Small cosmetic changes, removed 2nd derivative plot, included check 
% for "bad profiles" i.e. not flake profiles. 5/23/2018
% v1.0: First implementation of code. 5/16/2018
% David Lam, (c)2018

%% Data Import 
[filepath,name,ext] = fileparts(fileName);
fileID = fopen(fileName);
currentLine = fgetl(fileID);
data = [];
while currentLine ~= -1
    currentRow = textscan(currentLine, '%s');
    data = [data; str2double(currentRow{1,1})'];
    currentLine = fgetl(fileID);
end
fclose(fileID);

%% Profile Analysis
% Constants that are used in the rest of the code. These should not be
% changed unless you are aware of how they work. The numbers are also pretty
% arbitrary and are emperically determined (aka I guessed). 

% numProfiles: the number of profiles in a given imported file
% pointsInterp: number of points that goes into the smoothing of an
%               indivdual profile. Default is 1,000
% pointsRange: number of points to look past the curve maximum both left and
%              right for a cutoff in finding the minimum second derivative.
%              Default is 20% of the interpolated points
% pointsRange2ndDeriv: number of points to look past the 2nd derivative
%                      maxima both left and right in finding the appropriate
%                      cutoff point. Default is 25% of the interpolated
%                      points. 
% pointsRangeEdge: number of points to look away from the edges of the
%                  profile if the 2nd derivative maxima is at the edge.
%                  Default is 2% of the interpolated points. 
% cutoff: the value where you determine flake edges as a function of the
%         value of the maximum second derivative on each side. Default is
%         1/e. 
numProfiles = length(data(1,:))/2;
pointsInterp = 1000; 
pointsRange= 0.20*pointsInterp; 
pointsRange2ndDeriv = 0.25*pointsInterp; 
pointsRangeEdge = 0.02*pointsInterp; 
cutoff = exp(-1); 

% Initializes the matrices for thickness and lengths for histograms. 
thicknesses = zeros(numProfiles,1);
lengths = zeros(numProfiles,1);

% Begins data analysis on each profile. 
if(programMode ~= "area")
    for prof = 1:numProfiles
        % Boolean to include the final data. Can be modified via programMode.
        includeData = 1; 

        % Data sanitization to remove NaNs
        x = data(:, 2*prof-1);
        x = x(~isnan(x));
        y = data(:, 2*prof);
        y = y(~isnan(y));

        % Smoothing data out with spline fitting
        xinterp = linspace(x(1), x(end), pointsInterp);
        curve = interp1(x,y,xinterp, 'spline');
        [curveMax, curveMaxPos] = max(curve);

        % Take second derivative to find max negative second derivative and
        % onset of inflection points within some range about the curve maximum.
        secondD = diff(diff(curve));
        if curveMaxPos-pointsRange <= 0 || curveMaxPos+pointsRange > length(curve)
            warning('Check profile number %.0f for concavity issues!', prof)
            continue
        end
        minSecondD = min(secondD(curveMaxPos-pointsRange:curveMaxPos+pointsRange));
        minSecondDPos = find(secondD == minSecondD);

        % Divide the second derivative into a portion to the left and to the 
        % right of the determined second derivative min, then find the maximum. 
        xinterpSecondD = xinterp(2:end-1);
        leftSecondD = secondD(1:minSecondDPos);
        rightSecondD = secondD(minSecondDPos:end);
        [leftMax, indL] = max(leftSecondD);
        [rightMax, indR] = max(rightSecondD);

        % Looks at a subset of the left and right portion that range by the
        % constant pointsRange2ndDeriv
        leftSecondD = leftSecondD(max(1, indL-pointsRange2ndDeriv):indL);
        rightSecondD = rightSecondD(indR:min(indR+pointsRange2ndDeriv,length(rightSecondD)));

        % Finds the cutoff point for calculating flake length. 
        try
            leftPoint = find(secondD == max(leftSecondD(leftSecondD < cutoff*leftMax)));
        catch
            leftSecondD = secondD(1:pointsRangeEdge); 
            leftPoint = find(secondD == min(leftSecondD));
            warning('Check profile number %.0f for possible left edge issues!', prof);
            if(programMode == "noErrors") 
                includeData = 0; 
            end
        end

        try
            rightPoint = find(secondD == max(rightSecondD(rightSecondD < cutoff*rightMax)));
        catch
            rightSecondD = secondD(length(secondD)-pointsRangeEdge:length(secondD));
            rightPoint = find(secondD == min(rightSecondD));
            warning('Check profile number %.0f for possible right edge issues!', prof);
            if(programMode == "noErrors") 
                includeData = 0; 
            end
        end

        if(programMode == "manual")
            % Plotting the profile data and relevant bounds. 
            figure('units','normalized','outerposition',[0 0 1 1])
            yyaxis left
            plot(xinterp, curve, 'LineWidth', 3)
            hold on
            title(strcat(name, " Profile Number ", num2str(prof)), 'Interpreter', 'none')
    %         plot(xinterp(curveMaxPos), curveMax, 'g.', 'MarkerSize', 24)
    %         plot([xinterp(curveMaxPos-pointsRange), xinterp(curveMaxPos-pointsRange)], ylim, '-.', 'LineWidth', 3)
    %         plot([xinterp(curveMaxPos+pointsRange), xinterp(curveMaxPos+pointsRange)], ylim, '-.', 'LineWidth', 3)
            hold off

            % Plotting the second derivative. 
            yyaxis right
            hold on 
    %         plot(xinterpSecondD, secondD, 'LineWidth', 3);
    %         plot([xinterpSecondD(1), xinterpSecondD(end)], [0, 0],'-')
    %         plot(xinterpSecondD(minSecondDPos), secondD(minSecondDPos), 'k.', 'MarkerSize', 24)
    %         plot(xinterpSecondD(indL), secondD(indL), 'k.', 'MarkerSize', 24)
    %         plot(xinterpSecondD(indR+minSecondDPos), secondD(indR+minSecondDPos), 'k.', 'MarkerSize', 24)
            plot([xinterpSecondD(leftPoint), xinterpSecondD(leftPoint)], ylim, '-.', 'LineWidth', 3)
            plot([xinterpSecondD(rightPoint), xinterpSecondD(rightPoint)], ylim, '-.', 'LineWidth', 3)

            hold off

            % Prompts user if they want to include the profile. If not,
            % thickness and length with not be included in final output. 
            userAnswer = MFquestdlg([1,0.5],'Include this profile in your analysis?', 'Manual Profile Analysis');
            close
            if userAnswer == "No"
                includeData = 0;
            elseif userAnswer == "Cancel"
                break
            end
        end
        
        % Only adds data if the data is deemed to be included. "auto" includes
        % all data, while "manual" mode lets you pick the profiles for inclusion
        % and "noErrors" removes anything with edge errors. 
        if(includeData)
            % Calculates the thickness and length of a flake from the given profile in nm.
            thickness = (max(curve(leftPoint:rightPoint))-min(curve(leftPoint:rightPoint)))*10^9;
            len = (xinterpSecondD(rightPoint) - xinterpSecondD(leftPoint))*10^9;

            % Includes calculated values into the data. 
            thicknesses(prof) = thickness;
            lengths(prof) = len;
        end
    end
else
    for pair = 1:numProfiles/2
        % Boolean to include the final data. Can be modified via programMode.
        includeData = 1; 
        % Initializes a pair of thickness and length
        thickness = [0 0];
        len = [0 0];
        for i = 1:2
            prof = 2*pair -2 + i; %Selects the first, then second profile in a pair. 
            x = data(:, 2*prof-1);
            x = x(~isnan(x));
            y = data(:, 2*prof);
            y = y(~isnan(y));

            % Smoothing data out with spline fitting
            xinterp = linspace(x(1), x(end), pointsInterp);
            curve = interp1(x,y,xinterp, 'spline');
            [curveMax, curveMaxPos] = max(curve);

            % Take second derivative to find max negative second derivative and
            % onset of inflection points within some range about the curve maximum.
            secondD = diff(diff(curve));
            if curveMaxPos-pointsRange <= 0 || curveMaxPos+pointsRange > length(curve)
                warning('Check profile number %.0f for concavity issues!', prof)
                continue
            end
            minSecondD = min(secondD(curveMaxPos-pointsRange:curveMaxPos+pointsRange));
            minSecondDPos = find(secondD == minSecondD);

            % Divide the second derivative into a portion to the left and to the 
            % right of the determined second derivative min, then find the maximum. 
            xinterpSecondD = xinterp(2:end-1);
            leftSecondD = secondD(1:minSecondDPos);
            rightSecondD = secondD(minSecondDPos:end);
            [leftMax, indL] = max(leftSecondD);
            [rightMax, indR] = max(rightSecondD);

            % Looks at a subset of the left and right portion that range by the
            % constant pointsRange2ndDeriv
            leftSecondD = leftSecondD(max(1, indL-pointsRange2ndDeriv):indL);
            rightSecondD = rightSecondD(indR:min(indR+pointsRange2ndDeriv,length(rightSecondD)));

            % Finds the cutoff point for calculating flake length. 
            try
                leftPoint = find(secondD == max(leftSecondD(leftSecondD < cutoff*leftMax)));
            catch
                leftSecondD = secondD(1:pointsRangeEdge); 
                leftPoint = find(secondD == min(leftSecondD));
                warning('Check profile number %.0f for possible left edge issues!', prof);
            end

            try
                rightPoint = find(secondD == max(rightSecondD(rightSecondD < cutoff*rightMax)));
            catch
                rightSecondD = secondD(length(secondD)-pointsRangeEdge:length(secondD));
                rightPoint = find(secondD == min(rightSecondD));
                warning('Check profile number %.0f for possible right edge issues!', prof);
            end

            % Plotting the profile data and relevant bounds. 
            figure('units','normalized','outerposition',[0 0 1 1])
            yyaxis left
            plot(xinterp, curve, 'LineWidth', 3)
            hold on
            title(strcat(name, " Profile Number ", num2str(prof)), 'Interpreter', 'none')
    %         plot(xinterp(curveMaxPos), curveMax, 'g.', 'MarkerSize', 24)
    %         plot([xinterp(curveMaxPos-pointsRange), xinterp(curveMaxPos-pointsRange)], ylim, '-.', 'LineWidth', 3)
    %         plot([xinterp(curveMaxPos+pointsRange), xinterp(curveMaxPos+pointsRange)], ylim, '-.', 'LineWidth', 3)
            hold off

            % Plotting the second derivative. 
            yyaxis right
            hold on 
    %         plot(xinterpSecondD, secondD, 'LineWidth', 3);
    %         plot([xinterpSecondD(1), xinterpSecondD(end)], [0, 0],'-')
    %         plot(xinterpSecondD(minSecondDPos), secondD(minSecondDPos), 'k.', 'MarkerSize', 24)
    %         plot(xinterpSecondD(indL), secondD(indL), 'k.', 'MarkerSize', 24)
    %         plot(xinterpSecondD(indR+minSecondDPos), secondD(indR+minSecondDPos), 'k.', 'MarkerSize', 24)
            plot([xinterpSecondD(leftPoint), xinterpSecondD(leftPoint)], ylim, '-.', 'LineWidth', 3)
            plot([xinterpSecondD(rightPoint), xinterpSecondD(rightPoint)], ylim, '-.', 'LineWidth', 3)

            hold off

            % Prompts user if they want to include the profile. If not,
            % thickness and length with not be included in final output. 
            userAnswer = MFquestdlg([0.8,0.5],'Include this profile in your analysis?', 'Manual Profile Analysis');
            close
            if userAnswer == "No"
                includeData = 0;
                break;
            elseif userAnswer == "Cancel"
                return
            end
            thickness(i) = (max(curve(leftPoint:rightPoint))-min(curve(leftPoint:rightPoint)))*10^9;
            len(i) = (xinterpSecondD(rightPoint) - xinterpSecondD(leftPoint))*10^9;
        end
        
        % Only adds data if the pair of data is deemed to be included. "auto" includes
        % all data, while "manual" mode lets you pick the profiles for inclusion
        % and "noErrors" removes anything with edge errors.
        userAnswer = MFquestdlg([1,0.5],sprintf('Are your two profiles representative of the same flake?\nThickness 1: %3.2f nm, Thickness 2: %3.2f nm, Delta T: %3.2f nm'...
            , [thickness(1), thickness(2), abs(thickness(1)-thickness(2))]), 'Confirm Flake Profiles?');
        drawnow; %Stops UI lag from the questdlg
        close
        if userAnswer == "Yes"
            thicknesses(2*pair-1:2*pair) = thickness;
            lengths(2*pair -1:2*pair) = len;
        elseif userAnswer == "No"
            continue
        elseif userAnswer == "Cancel"
            break
        end
    end
end



% Removes non-zero values for the thicknesses, though this can be
% surpressed to find the profile numbers that have "failed". 
notIncludedIndexes = find(thicknesses == 0)
thicknesses = thicknesses(thicknesses ~= 0); 
lengths = lengths(lengths ~=0); 

if(programMode ~= "area") 
    vertDim = thicknesses;
    latDim = lengths; 
end

if(programMode == "area")
    vertDim = zeros(size(thicknesses,1)/2, 1);
    latDim = zeros(size(lengths,1)/2, 1);
    for i = 1:size(thicknesses,1)/2
        p1 = 2*i-1;
        p2 = 2*i;
        vertDim(i) = mean(thicknesses(p1:p2));
        latDim(i) = areaFunction(lengths(p1),lengths(p2)); %Invokes some model of flake area given 2 distances.
    end
end
end

function [flakeArea] = areaFunction(a1, a2, model)
    flakeArea = pi*a1*a2/4; %Elliptical model 
end

