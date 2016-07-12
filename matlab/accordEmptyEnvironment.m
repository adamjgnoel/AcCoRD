function [hFig, hAxes] = accordEmptyEnvironment(fileToLoad, scale, ...
    customFigProp, customAxesProp, ...
    regionToPlot, customRegionProp, actorToPlot, customActorProp,...
    cameraAnchor)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordEmptyEnvironment.m - plot the regions and actors in an AcCoRD
%   configuration file. Not intended for drawing molecules (use
%   accordVideoMaker to draw simulation observations onto environment)
%
% INPUTS
% fileToLoad - simulation file generated by accordImport
% scale - scaling of physical dimensions of region and actor coordinates.
%   Needed to mitigate patch display problems. Recommend that smallest
%   object (non-molecule) to plot has dimension of order 1
% customFigProp - structure of figure properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildFigureStruct for structure fields and their default values.
% customAxesProp - structure of axes properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildAxesStruct for structure fields and their default values.
% regionToPlot - array of indices of regions to be plotted.
% customRegionProp - structure of region properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildDispStruct for structure fields and their default values.
% actorToPlot - array of indices of actors to be plotted. Indexing matches
%   the actor list in the original config file and is independent of
%   whether an actor is active or passive. Actors listed here will have
%   their shapes plotted but NOT their molecules (the latter is indicated
%   by the argument "passiveActorToPlot". Actors are drawn after regions,
%   so if an actor is defined by region(s) then it will be drawn on top of
%   its region(s).
% customActorProp - structure of actor properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildDispStruct for structure fields and their default values.
% cameraAnchor - cell array defining camera display anchor.
%   Can be passed as an empty cell array. Defines camera settings, in
%   the format {'CameraPosition', 'CameraTarget', 'CameraViewAngle',
%   'CameraUpVector'}. See MATLAB camera documentation for more details.
%
% OUTPUTS
% hFig - handle(s) to plotted figure(s). Use for making changes.
% hAxes - handle(s) to axes in plotted figure(s). Use for making changes.
%
% Last revised for AcCoRD v0.7 (public beta, 2016-07-09)
%
% Revision history:
%
% Revision v0.7 (public beta, 2016-07-09)
% - Created file
%
% Created 2016-06-28

%% Load Configuration File
% Config file is in JSON format
opt.SimplifyCell = 1;
opt.FastArrayParser = 1;
opt.ShowProgress = 0;
config = accordConfigImport(loadjson(fileToLoad, opt));

%% Load Default Display Properties and Apply Specified Changes
figureProp = accordBuildFigureStruct(customFigProp);
axesProp = accordBuildAxesStruct(customAxesProp);
regionDispStruct = accordBuildDispStruct(regionToPlot, customRegionProp);
actorDispStruct = accordBuildDispStruct(actorToPlot, customActorProp);

[hFig, hAxes] = accordPlotEnvironment(config,...
    axesProp, figureProp, regionDispStruct, actorDispStruct, scale);

if ~isempty(cameraAnchor)
    set(hAxes, {'CameraPosition','CameraTarget',...
        'CameraViewAngle','CameraUpVector'}, cameraAnchor);
end