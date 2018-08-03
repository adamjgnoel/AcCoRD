function [hFig, hAxes] = accordQuickPlot(filename)
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
%   the user modification required to quickly view results. For more
%   control, copy and modify accordPlotMakerWrapper.m
%
% INPUTS
% filename - .mat-file of AcCoRD simulation output, as saved by a call to
%   accordImport.m
%
% OUTPUTS
% hFig - handle to plotted figure. Use for making changes.
% hAxes - handle to axes in plotted figure. Use for making changes.
%
% Last revised for AcCoRD LATEST_VERSION
%
% Revision history:
%
% Revision LATEST_VERSION
% - Created file
%
% Created 2018-08-03

% Load AcCoRD Output
load(filename);

customFigProp = [];
% Change axes properties to favor curve plotting over video
customAxesProp.Visible = 'on';
customAxesProp.Projection = 'orthographic';
customAxesProp.Clipping = 'on';

customObsProp = [];

hAxes = 0;
numLines = 1;
colors = {'k', 'b', 'r', 'g', 'm', 'c'};
for i = 1:data.numPassiveRecord
    customCurveProp.DisplayName = ['Actor ' num2str(data.passiveRecordID(i))];
    for j = 1:data.passiveRecordNumMolType(i)
        customCurveProp.Color = colors{mod(numLines,length(colors))};
        [hFig, hAxes] = accordPlotMaker(hAxes, filename,...
            i, j, ...
            customObsProp, customCurveProp, ...
            customFigProp, customAxesProp);
        numLines = numLines + 1;
    end
end