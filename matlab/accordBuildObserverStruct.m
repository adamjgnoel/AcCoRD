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
%   what data to plot from a passive actor. Used by accordPlotMaker
%
% INPUTS
% propChange - cell array of structure of properties to modify from their
%   default values. Any field should match one in obsSpec
%
% OUTPUTS
% obsSpec - structure with observer display parameters
%
% %% HOWTO - Types of Plots Available %%
%   Change type of plot by assigning string to propChange.obsType. Options
%   and how to use them are below. All plots that use simulation data will
%   sample the observations stored in the simulation output file from the
%   indices 'firstSample':'sampleInterval':'lastSample', where 'lastSample' is
%   replaced with 'end' if 'bMaxSample' == true. The simulation realizations
%   to use are set via 'iterType' and (optionally) 'iterCustom'.
%   All plot types that do not use simulation data will ignore the sampling
%   and realization properties
% 
% 'Sample'
% - 2D curve of (average) time-varying simulation observations
% - set 'iterType' to 'Custom' and 'iterCustom' to a single index to
% plot a single realization (i.e., without averaging)
%
% 'Expected', '3D Expected'
% - 2D curve or 3D surface of custom data (non-simulation)
% - set 'data1' and 'data2' as vectors for 2D
% - set 'data1', 'data2', and 'data3' as matrices for 3D
%
% 'Empirical CDF', '3D Empirical CDF'
% - Empirical cumulative distribution function from simulation observations
% - no other unique structure members need to be set
% - 2D version will include all sampling times in one CDF (i.e., good for
% either 1 time only or for a system with constant average behavior)
% - 3D version will plot a separate CDF for each sampling time
%
% 'Expected CDF', '3D Expected CDF'
% - plot cumulative distribution function for specified probability
% distribution (non-simulation)
% - set distribution to 'Binomial', 'Poisson', or 'Gaussian'. Latter 2 are
% used as approximations of Binomial
% - set 'numTrials' and 'trialProbability' to define the number of trials
% associated with a given observation and its probability of success. In
% context of AcCoRD, number of molecules is the number of trials and
% probability of success is probability of given molecule being observed.
% - 2D version needs 'trialProbability' to be scalar
% - 3D version needs 'trialProbability' to be a vector, and 'data2' to be the
% same length where it lists the event times associated with the
% corresponding probabilities (i.e., we assume a system with uniform
% time-varying success probability). This version also uses data1 to define
% the CDF "bin" values (i.e., the measurement points
% in terms of # of successful trials). If not defined, range of bin values
% is generated from expected mean and variance
%
% 'Empirical PMF', '3D Empirical PMF'
% - Empirical probability mass function from simulation observations
% - no other unique structure members need to be set
% - 2D version will include all sampling times in one PMF (i.e., good for
% either 1 time only or for a system with constant average behavior)
% - 3D version will plot a separate PMF for each sampling time
%
% 'Expected PMF', '3D Expected PMF'
% - plot probability mass function for specified probability
% distribution (non-simulation)
% - set 'distribution' to 'Binomial', 'Poisson', or 'Gaussian'. Latter 2 are
% used as approximations of Binomial
% - set 'numTrials' and 'trialProbability' to define the number of trials
% associated with a given observation and its probability of success. In
% context of AcCoRD, number of molecules is the number of trials and
% probability of success is probability of given molecule being observed.
% - set data1 to define the PMF "bin" values (i.e., the measurement points
% in terms of # of successful trials). If not defined, range of bin values
% is generated from expected mean and variance
% - 2D version needs 'trialProbability' to be scalar
% - 3D version needs 'trialProbability' to be a vector, and 'data2' to be the
% same length where it lists the event times associated with the
% corresponding probabilities (i.e., we assume a system with uniform
% time-varying success probability)
%
% 'Mutual Information', '3D Mutual Information'
% - plot mutual information of simulation samples relative to 1 (2D) or
% more (3D) reference simulation samples
% - in 2D, first sample index is the "reference" and all mutual
% information calculations are relative to the first sample index. So, the
% first point on the curve will be the "perfect" case and measure the
% entropy of the first sample
% - in 3D, all sampled observations are "reference" observations. Set
% 'data1' as a vector to define the observations relative to each sample to
% calculate the mutual information between.
% The relative values here in 'data1' are applied to the original
% simulation data and are NOT scaled by the 'sampleInterval' parameter (the
% idea being that the scale of the offsets considered can be much smaller
% than the time between reference observations). If a value in 'data1' is
% 0, then the corresponding surface points will measure the entropy of each
% reference index.
%
% 'Monte Carlo Mutual Information', '3D Monte Carlo Mutual Information'
% - measure the "mutual information" between independent random variables
% (non-simulation). Theoretically this measurement should always be 0 but it
% is non-zero in practice for a finite number of samples. It is used to
% compare with the simulation mutual information to tell when the
% observations are "effectively indepedent".
% - all random variables are generated using Binomial distribution
% - set 'numTrials' and 'trialProbability' to define the number of trials
% associated with a given observation and its probability of success. In
% context of AcCoRD, number of molecules is the number of trials and
% probability of success is probability of given molecule being observed.
% - set data1 as vector to list x-axis values (would be reference times in
% corresponding simulation curve). In 2D, this can be any number of values
% if 'trialProbability' has length less than 3 since only 1 calculation will
% be made and the same value will be plotted for each point.
% - in 2D, set 'data2' to format [val1, val2], where val1 is the number of
% random numbers to generate (i.e., observe) for each variable in the
% mutual information calculation. val2 is the number of times to repeat the
% calculation (plotted result will be average of these calculations).
% - also in 2D, 'trialProbability' must be a vector where the first value
% is the "reference" observation probability. All mutual information
% calculations will be relative to random variables generated from this probability.
% - also in 2D, if there are exactly 2 values in 'trialProbability', then
% there will be only 1 mutual information determined (i.e., entropy by
% comparing first sample with itself is ignored).
% - in 3D, set 'data3' to format [val1, val2], where val1 is the number of
% random numbers to generate (i.e., observe) for each variable in the
% mutual information calculation. val2 is the number of times to repeat the
% calculation (plotted result will be average of these calculations).
% - also in 3D, 'trialProbability' must be a matrix. In each row, the first
% probability will be the "reference" observation probability that the
% remainder of the row is compared with. Set 'data2' to be the vector of
% y-axis values. 'trialProbability' must have length('data1') rows and
% length('data2') columns.
%
% Last revised for AcCoRD v0.7 (public beta, 2016-07-09)
%
% Revision history:
%
% Revision v0.7 (public beta, 2016-07-09)
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
    'iterType', 'All', ...                          % Include what simulation iterations. Options are 'All' and 'Custom'
    'iterCustom', 1);                               % Values to use if iterType == 'Custom'

%% Make Specified Changes to Defaults
if ~isempty(propChange)
    propFields = fieldnames(propChange);
    numProp = numel(propFields);
    for j = 1:numProp
        obsSpec.(propFields{j}) = propChange.(propFields{j});
    end
end