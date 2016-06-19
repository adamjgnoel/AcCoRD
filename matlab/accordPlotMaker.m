function [hFig, hAxes, hCurve] = accordPlotMaker(hAxes, fileToLoad,...
    passiveID, molID, customObsProp, customCurveProp, customFigProp, customAxesProp)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordPlotMaker.m - 
%
% INPUTS
%
% OUTPUTS
% hFig - handle(s) to plotted figure(s). Use for making changes.
% hAxes - handle(s) to axes in plotted figure(s). Use for making changes.
%
% Last revised for AcCoRD LATEST_VERSION
%
% Revision history:
%
% Revision LATEST_VERSION
% - Created file
%
% Created 2016-06-03


%% Load Default Display Properties and Apply Specified Changes
obsSpec = accordBuildObserverStruct(customObsProp);
curveSpec = accordBuildCurveStruct(customCurveProp);

%% Create figure and apply properties if it does not exist
if hAxes == 0
    figureProp = accordBuildFigureStruct(customFigProp);
    axesProp = accordBuildAxesStruct(customAxesProp);
    
    hFig = figure;
    if ~isempty(figureProp)
        figurePropFields = fieldnames(figureProp);
        numProp = numel(figurePropFields);
        for i = 1:numProp
            set(hFig, figurePropFields{i}, figureProp.(figurePropFields{i}));
        end
    end

    %% Create axes and apply properties
    hAxes = axes('Parent',hFig);
    if ~isempty(axesProp)
        axesPropFields = fieldnames(axesProp);
        numProp = numel(axesPropFields);
        for i = 1:numProp
            set(hAxes, axesPropFields{i}, axesProp.(axesPropFields{i}));
        end
    end
    hold(hAxes, 'on');
    
    switch obsSpec.obsType
        case 'Sample'
            xlabel(hAxes, 'Time [s]');
            ylabel(hAxes, 'Number of molecules');
        case 'Empirical CDF'
            xlabel(hAxes, 'Number of Molecules');
            ylabel(hAxes, 'Cumulative Distribution Function');
        case 'Histogram'
            xlabel(hAxes, 'Number of Molecules');
            ylabel(hAxes, 'Frequency');
    end
else
    hFig = get(hAxes, 'Parent');
end


%% Add Curve to plot
load(fileToLoad, 'config', 'data');

obsMatrix = data.passiveRecordCount{passiveID}(:,molID,:);
numObs = size(obsMatrix,3);
numRepeat = size(obsMatrix,1);

% Determine all original simulation observation times
actorID = 1;
while config.actor{actorID}.passiveID ~= passiveID
    actorID = actorID + 1;
end
tArrayFull = config.actor{actorID}.startTime + ...
    (0:(numObs-1))*config.actor{actorID}.actionInterval;

% Determine times that will be sampled
if obsSpec.bMaxSample
    tPlotInd = obsSpec.firstSample:obsSpec.sampleInterval:numObs;
else
    tPlotInd = obsSpec.firstSample:obsSpec.sampleInterval:obsSpec.lastSample;
end
obsMatrixSampled = obsMatrix(:,:,tPlotInd);

% Determine data to plot
switch obsSpec.obsType
    case 'Sample'
        obsPlotMatrix = obsMatrixSampled;
        xData = tArrayFull(tPlotInd);
        % Determine what repeats to average over (default is all)
        switch obsSpec.avgType
            case 'All'
                avgInd = 1:numRepeat;
            case 'Custom'
                avgInd = obsSpec.avgCustom;
            otherwise
        end
        yData = reshape(mean(obsPlotMatrix(avgInd,:,:),1),1,[]);
    case 'Empirical CDF'
        [yData, xData] = ecdf(obsMatrixSampled(:));
    case 'Histogram'
        % yData and xData will be calculated later
        [yData, xData] = histcounts(obsMatrixSampled(:), obsSpec.numHistBins);
        numX = length(xData);
        xData = (xData(1:(numX-1))+xData(2:numX))/2;
        i = 1;
        while i < numX
            if yData(i) == 0
                yData(i) = [];
                xData(i) = [];
                numX = numX - 1;
            else
                i = i + 1;
            end
        end
        
            
    otherwise
        error('Observation type "%s" invalid', obsType)
end

% Normalize the x-data?
if obsSpec.bNormalizeX
    % Normalize with respect to what?
    switch obsSpec.normalizeTypeX
        case 'Max'
            maxValX = max(yData);
        case 'Custom'
            maxValX = obsSpec.normalizeCustomX;
        otherwise
    end
else
    maxValX = 1;
end

% Normalize the y-data?
if obsSpec.bNormalizeY
    % Normalize with respect to what?
    switch obsSpec.normalizeTypeY
        case 'Max'
            maxValY = max(yData);
        case 'Custom'
            maxValY = obsSpec.normalizeCustomY;
        otherwise
    end
else
    maxValY = 1;
end

%% Plot Curve and Format it
hCurve = plot(hAxes, xData./maxValX, yData./maxValY);
if ~isempty(curveSpec)
    curvePropFields = fieldnames(curveSpec);
    numProp = numel(curvePropFields);
    for i = 1:numProp
        set(hCurve, curvePropFields{i}, curveSpec.(curvePropFields{i}));
    end
end

%% Cleanup
end