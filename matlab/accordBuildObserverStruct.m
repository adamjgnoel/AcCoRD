function obsSpec = accordBuildObserverStruct(propChange)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordBuildObserverStruct.m - build structure with parameters to define
%   what data to plot from a passive actor
%
% INPUTS
% propChange - cell array of structure of properties to modify from their
%   default values. Any field must match one in obsSpec and be the
%   same size.
%
% OUTPUTS
% obsSpec - structure with observer display parameters
%
% Last revised for AcCoRD LATEST_VERSION
%
% Revision history:
%
% Revision LATEST_VERSION
% - Created file
%
% Created 2016-06-03

%% Set Default Values
obsSpec = struct('obsType', 'Sample', ...           % Type of data to plot. Options are 'Sample',
    ...                                             % 'Empirical CDF', '3D Empirical CDF',
    ...                                             % 'Histogram', and '3D Histogram'
    'numHistBins', 100, ...                         % # of histogram bins if obsType == 'Histogram' or '3D Histogram'
    'firstSample', 1, ...                           % Index of first observation to sample
    'sampleInterval', 1, ...                        % Frequency of observation sampling
    'bMaxSample', true, ...                         % Keep sampling until last observation?
    'lastSample', 1, ...                            % Index of last observation to sample (if bMaxSample == false)
    'bNormalizeX', false, ...                       % Normalize x-data?
    'normalizeTypeX', 'Max', ...                    % Type of x-normalization. Options are 'Max' and 'Custom'
    'normalizeCustomX', 1, ...                      % Value to use if normalizeTypeX == 'Custom'
    'bNormalizeY', false, ...                       % Normalize y-data?
    'normalizeTypeY', 'Max', ...                    % Type of y-normalization. Options are 'Max' and 'Custom'
    'normalizeCustomY', 1, ...                      % Value to use if normalizeTypeY == 'Custom'
    'bNormalizeZ', false, ...                       % Normalize z-data? Only applies for 3D plots
    'normalizeTypeZ', 'Max', ...                    % Type of z-normalization. Options are 'Max' and 'Custom'
    'normalizeCustomZ', 1, ...                      % Value to use if normalizeTypeZ == 'Custom'
    'avgType', 'All', ...                           % Average results over what iterations. Options are 'All' and 'Custom'
    'avgCustom', 1);                                % Values to use if avgType == 'Custom'

%% Make Specified Changes to Defaults
if ~isempty(propChange)
    propFields = fieldnames(propChange);
    numProp = numel(propFields);
    for j = 1:numProp
        obsSpec.(propFields{j}) = propChange.(propFields{j});
    end
end