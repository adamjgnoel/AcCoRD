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
%   default values. Any field should match one in obsSpec
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
    ...                                             %   'Empirical CDF', '3D Empirical CDF',
    ...                                             %   'Empirical PMF', '3D Empirical PMF',
    ...                                             %   'Expected', 'Expected PMF', '3D Expected PMF',
    ...                                             %   'Expected CDF', '3D Expected CDF',
    ...                                             %   'Mutual Information', '3D Mutual Information'
    ...                                             %   'Monte Carlo Mutual Information', '3D Monte Carlo Mutual Information'
    'distribution', 'Binomial', ...                 % Statistical distribution if obsType is an Expected CDF or PMF
    ...                                             %   or expected CDF
    'trialProbability', 0, ...                      % Event probability for Expected CDF or PDF. Should be a vector for 3D versions
    'numTrials', 0, ...                             % Number of event trails for Expected CDF or PDF
    ...                                             %   or Monte Carlo Mutual Information
    'numMeasurements', 100, ...                     % # of bins if obsType == 'Expected CDF', '3D Expected CDF'
    ...                                             %   # of repeated mutual information measurements for Monte Carlo Mutual Information
    'data1', [], ...                                % Custom data vector 1
    ...                                             %   Defines xData if obsType == 'Expected'
    ...                                             %   Overwrites distribution bin values if obsType is an expected PMF or CDF
    ...                                             %   Determines plotting x-coordinate if obsType == 'Monte Carlo Mutual Information'
    'data2', [], ...                                % Custom data vector 2
    ...                                             %   Defines yData if obsType == 'Expected'
    'data3', 1, ...                                 % Custom data vector 3
    ...                                             %   Also used for superposition of 2D data over 3D surfaces
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