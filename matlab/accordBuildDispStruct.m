function dispStruct = accordBuildDispStruct(indToDisp)
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
%   some kind of object should be plotted (i.e., Region or actor)
%
% INPUTS
% indToDisp - indices of objects that will be displayed and thus need
%   properties
%
% OUTPUTS
% dispStruct - structure with display parameters
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

dispStruct = struct('numToDisp', indLength, ...
    'indToDisp', indToDisp, ...
    'bDispFace', false(1, indLength), ...
    'dispColor', cell(1), ...
    'opaque', ones(1, indLength));

% Cell array needs to be defined as a cell with a cell array in order to
% suppress the output structure being an array of structures
dispStruct.dispColor{1} = cell(1,indLength);

for i = 1:indLength
    dispStruct.dispColor{1}{i} = 'black';
end