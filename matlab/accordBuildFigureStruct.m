function figureProp = accordBuildFigureStruct(xPixels, yPixels)
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
%   are considered "defaults".
%
% INPUTS
% xPixels - size of figure along x direction (in pixels). Optional (can be
%   empty)
% yPixels - size of figure along y direction (in pixels). Optional (can be
%   empty)
%
% OUTPUTS
% figureProp - structure with common figure properties
%
% Last revised for AcCoRD LATEST_VERSION
%
% Revision history:
%
% Revision LATEST_VERSION
% - Created file
%
% Created 2016-05-18

% Overall axes placement and visibility
figureProp.Color = 'w';     % If on, draw complete box outline along axis
                            % limits
figureProp.Clipping = 'on'; % If on, crop parts of objects that are defined
                            % beyond figure limits
figureProp.Position = [100 100 560 420]; % Relative position of figure on screen
if ~isempty(xPixels)
    figureProp.Position(3) = xPixels;
end
if ~isempty(yPixels)
    figureProp.Position(4) = yPixels;
end