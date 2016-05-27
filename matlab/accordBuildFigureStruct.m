function figureProp = accordBuildFigureStruct(propChange)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordBuildFigureStruct.m - build structure with desired parameters for
%   plotting figure. All structure fields must be valid figure properties
%   since they will be directly applied to the axes handle via set.
%   Field values can be appended or changed within this function or by
%   modifying the output structure. The fields and parameters listed here
%   are considered "defaults". By default, the figure is drawn white,
%   objects drawn beyond the figure limits are clipped, and the default
%   MATLAB figure size is used.
%
% INPUTS
% propChange - structure of properties to modify from their default values
%   Any field must be a standard Matlab figure property and be the correct
%   size. Fields that match any of the default values below will be applied
%   in the order below. Fields that do not match the defaults below will
%   be applied after the defaults and in the order defined.
%
% OUTPUTS
% figureProp - structure with common figure properties
%
% Last revised for AcCoRD v0.6 (public beta, 2016-05-30)
%
% Revision history:
%
% Revision v0.6 (public beta, 2016-05-30)
% - Created file
%
% Created 2016-05-18

%% Set Default Values
% Overall axes placement and visibility
figureProp.Color = 'w';     % If on, draw complete box outline along axis
                            % limits
figureProp.Clipping = 'on'; % If on, crop parts of objects that are defined
                            % beyond figure limits
figureProp.Position = [100 100 560 420]; % Relative position of figure on screen

%% Make Specified Changes to Defaults
if ~isempty(propChange)
    propFields = fieldnames(propChange);
    numProp = numel(propFields);
    for i = 1:numProp
        figureProp.(propFields{i}) = propChange.(propFields{i});
    end
end