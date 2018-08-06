function [hFig, hAxes] = accordQuickPlot(inputFile)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2018 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordQuickPlot.m - Wrapper function for accordPlotMaker, which
%   can create a plot from AcCoRD simulation output. Uses defaults to limit
%   the user modification that is required. Plots the (average)
%   time-varying of ALL types of molecules that were observed by ALL
%   passive actors. All curves are plotted as solid lines, and the colors
%   cycle through black, blue, red, green, magenta, and cyan
%   For more control and options, copy and modify accordPlotMakerWrapper.m.
%
% INPUTS
% filename - full filename (including directory) of data to plot. This can
%   be either the "root" filename of simulation output (i.e., the txt file
%   without the "_SEED##" suffix), which MUST have seed number 1, or a
%   '.mat'-file of output that was already imported via a call to
%   accordImport.m. Both 'txt' and 'mat' file extensions are optional.
%
% OUTPUTS
% hFig - handle to plotted figure. Use for making changes.
% hAxes - handle to axes in plotted figure. Use for making changes.
%
% Last revised for AcCoRD v1.4 (2018-08-06)
%
% Revision history:
%
% Revision v1.4 (2018-08-06)
% - Created file
%
% Created 2018-08-03

% Determine whether filename refers directly to simulation output or to a
% simulation that was previously loaded
[filepath, fileName, fileEnd] = fileparts(inputFile);
if strcmp(fileEnd, '.mat')
    % We already have simulation output that we can load directly
    load(inputFile, 'data');
    dataFileName = inputFile;
elseif isempty(fileEnd) && strcmp(fileName((end-3):end), '_out')
    % We already have simulation output that we can load directly
    dataFileName = fullfile(filepath, strcat(fileName, '.mat'));
    load(dataFileName, 'data');
else
    % Assume we just have simulation data that we need to import first
    if strcmp(fileEnd, '.txt')
        inputFile = fullfile(filepath, fileName);
    end
    [data, ~] = accordImport(inputFile, 1, 1);
    dataFileName = strcat(fileName, '_out.mat');
end

customFigProp = [];
% Change axes properties to favor curve plotting over video
customAxesProp.Visible = 'on';
customAxesProp.Projection = 'orthographic';
customAxesProp.Clipping = 'on';

customObsProp = [];

hAxes = 0;
numLines = 1;
colors = {'k', 'b', 'r', 'g', 'm', 'c', 'g', 'y'};
for i = 1:data.numPassiveRecord
    for j = 1:data.passiveRecordNumMolType(i)
        customCurveProp.DisplayName = ['Actor ' num2str(data.passiveRecordID(i)) ...
            ' Molecule Type ' num2str(data.passiveRecordMolID{i}(j))];
        customCurveProp.Color = colors{mod(numLines-1,length(colors))+1};
        [hFig, hAxes] = accordPlotMaker(hAxes, dataFileName,...
            i, j, ...
            customObsProp, customCurveProp, ...
            customFigProp, customAxesProp);
        numLines = numLines + 1;
    end
end
legend(hAxes, 'show')