function markerStruct = accordBuildMarkerStruct(indToDisp)
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
% indToDisp - indices of molecule types that will be displayed and thus need
%   properties
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

indLength = length(indToDisp);

markerStruct = struct('numToDisp', indLength, ...   % # of types being plotted
    'indToDisp', indToDisp, ...                     % indices of those being plotted
    'dispStr', cell(1), ...                         % Display string (if legend turned on)
    'edgeColor', cell(1), ...                       % Marker edge (border) color
    'faceColor', cell(1), ...                       % Marker face color
    'marker', cell(1), ...                          % Marker shape
    'lineWidth', 0.5*ones(1, indLength), ...        % Marker edge (border) width
    'size', 50*ones(1, indLength));                 % Marker size in points^2, where 1 point = 1/72 inch

% Cell arrays need to be defined as a cells with cell arrays in order to
% suppress the output structure being an array of structures
markerStruct.dispStr{1} = cell(1,indLength);
markerStruct.edgeColor{1} = cell(1,indLength);
markerStruct.faceColor{1} = cell(1,indLength);
markerStruct.marker{1} = cell(1,indLength);

for i = 1:indLength
    markerStruct.dispStr{1}{i} = sprintf('Molecule %d',indToDisp(i));
    markerStruct.edgeColor{1}{i} = 'black';
    markerStruct.faceColor{1}{i} = 'none';
    markerStruct.marker{1}{i} = 'o';
end