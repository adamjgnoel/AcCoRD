function markerStruct = accordBuildMarkerStruct(numActor, indToDisp, propChange)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordBuildMarkerStruct.m - build structure with parameters to define how
%   markers should be plotted (i.e., molecule types)
%
% INPUTS
% numActor - number of actors whose observations are being recorded
% indToDisp - cell array of indices of molecule types that will be
%   displayed and thus need properties
% propChange - cell array of structure of properties to modify from their
%   default values. Any field must match one in markerStruct and be the
%   same size.
%
% OUTPUTS
% markerStruct - structure with marker display parameters
%
% Last revised for AcCoRD LATEST_VERSION
%
% Revision history:
%
% Revision LATEST_VERSION
% - Created file
%
% Created 2016-05-18

markerStruct = cell(1,numActor);

for i = 1:numActor
    %% Set Default Values

    indLength = length(indToDisp{i});

    markerStruct{i} = struct('numToDisp', indLength, ...   % # of types being plotted
        'indToDisp', indToDisp{i}, ...                     % indices of those being plotted
        'dispStr', cell(1), ...                         % Display string (if legend turned on)
        'edgeColor', cell(1), ...                       % Marker edge (border) color
        'faceColor', cell(1), ...                       % Marker face color
        'marker', cell(1), ...                          % Marker shape
        'lineWidth', 0.5*ones(1, indLength), ...        % Marker edge (border) width
        'size', 50*ones(1, indLength));                 % Marker size in points^2, where 1 point = 1/72 inch

    % Cell arrays need to be defined as a cells with cell arrays in order to
    % suppress the output structure being an array of structures
    markerStruct{i}.dispStr{1} = cell(1,indLength);
    markerStruct{i}.edgeColor{1} = cell(1,indLength);
    markerStruct{i}.faceColor{1} = cell(1,indLength);
    markerStruct{i}.marker{1} = cell(1,indLength);

    for j = 1:indLength
        markerStruct{i}.dispStr{1}{j} = sprintf('Molecule %d',indToDisp{i}(j));
        markerStruct{i}.edgeColor{1}{j} = 'black';
        markerStruct{i}.faceColor{1}{j} = 'none';
        markerStruct{i}.marker{1}{j} = 'o';
    end

    %% Make Specified Changes to Defaults
    if ~isempty(propChange{i})
        propFields = fieldnames(propChange{i});
        numProp = numel(propFields);
        for j = 1:numProp
            markerStruct{i}.(propFields{j}) = propChange{i}.(propFields{j});
        end
    end
end