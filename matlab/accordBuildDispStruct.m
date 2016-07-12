function dispStruct = accordBuildDispStruct(indToDisp, propChange)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordBuildDispStruct.m - build structure with parameters to define how
%   some kind of object should be plotted (i.e., Region or actor). Default
%   settings will plot only the edges of the object (and not the faces) in
%   fully opaque black.
%
% INPUTS
% indToDisp - indices of objects that will be displayed and thus need
%   properties
% propChange - structure of properties to modify from their default values.
%   Any field must match one in dispStruct and be the same size.
%
% OUTPUTS
% dispStruct - structure with display parameters
%
% Last revised for AcCoRD v0.7 (public beta, 2016-07-09)
%
% Revision history:
%
% Revision v0.7 (public beta, 2016-07-09)
% - added option to display region subvolumes. Will only apply if object is
%   a region; will be ignored for actors.
%
% Revision v0.6 (public beta, 2016-05-30)
% - Created file
%
% Created 2016-05-18

%% Set Default Values

indLength = length(indToDisp);

dispStruct = struct('numToDisp', indLength, ...
    'indToDisp', indToDisp, ...
    'bDispFace', false(1, indLength), ...
    'dispColor', cell(1), ...
    'opaque', ones(1, indLength), ...
    'bDispSubvolumes', false(1,indLength));

% Cell array needs to be defined as a cell with a cell array in order to
% suppress the output structure being an array of structures
dispStruct.dispColor{1} = cell(1,indLength);

for i = 1:indLength
    dispStruct.dispColor{1}{i} = 'black';
end

%% Make Specified Changes to Defaults
if ~isempty(propChange)
    propFields = fieldnames(propChange);
    numProp = numel(propFields);
    for i = 1:numProp
        dispStruct.(propFields{i}) = propChange.(propFields{i});
    end
end