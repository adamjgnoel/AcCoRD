function curveSpec = accordBuildSurfStruct(propChange)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordBuildSurfStruct.m - build structure with parameters to define how
%   to plot a surface of simulation data
%
% INPUTS
% propChange - cell array of structure of properties to modify from their
%   default values. Any field must match a Matlab Chart Surface Property
%
% OUTPUTS
% curveSpec - structure with surface parameters
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
curveSpec = struct(...
    'LineStyle', '-', ...           % Edge Line style. Options: '-', '--', ':', '-.', 'none'
    'LineWidth', 0.5, ...           % Edge Line width
    'EdgeColor', 'flat', ...        % Edge color. Options: 'flat','interp','none','texturemap', RGB triplet, color string
    'FaceColor', 'flat', ...        % Face color. Options: 'flat','interp','none','texturemap', RGB triplet, color string
    'Marker', 'none', ...           % Marker shape. Options: 'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h', 'none'
    'MarkerSize', 6, ...            % Marker size
    'MarkerEdgeColor', 'auto', ...  % Marker edge color. Options: 'auto', 'none', 'RGB triplet', 'color string'
    'MarkerFaceColor', 'none', ...  % Marker edge color. Options: 'auto', 'none', 'RGB triplet', 'color string'
    'DisplayName', '');             % Display string (if legend turned on)

%% Make Specified Changes to Defaults
if ~isempty(propChange)
    propFields = fieldnames(propChange);
    numProp = numel(propFields);
    for j = 1:numProp
        curveSpec.(propFields{j}) = propChange.(propFields{j});
    end
end