function [hFig, hAxes] = accordEmptyEnvironmentQuick(filename)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016-2019 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordEmptyEnvironmentQuick.m - Wrapper function for
%   accordEmptyEnvironment, which creates a figure displaying the regions
%   and actors described by an AcCoRD configuration file.
%   This file prepares all of the input arguments needed for
%   accordEmptyEnvironment, so that ALL actors and regions are plotted with
%   default settings.
%
% INPUTS
% NONE - this file is intended to be a wrapper that is modified directly.
%   Descriptions of the input arguments for accordEmptyEnvironment are
%   described as they are initialized below.
%
% OUTPUTS
% hFig - handle to plotted figure. Use for making changes.
% hAxes - handle to axes in plotted figure. Use for making changes.
%
% Last revised for AcCoRD v1.4.1 (2019-04-16)
%
% Revision history:
%
% Revision v1.4.1 (2019-04-16)
% - Created file from accordEmptyEnvironmentWrapper
% - modified default list of actors and regions to plot so that EVERYTHING
% is plotted and not just the first region and first actor. This default
% list is more user friendly.
%
% Created 2018-09-12

%% 

% fileToLoad - .txt file (with extension) of environment configuration
%   data. Include relative directory structure if applicable. File must be
%   a valid AcCoRD JSON-formatted configuration file in order to be
%   imported correctly.
fileToLoad = filename;

% Load the file so that we can use the full list of regions and actors as
% the defaults to plot
addpath('JSONlab');
opt.SimplifyCell = 1;
opt.FastArrayParser = 1;
opt.ShowProgress = 0;
config = accordConfigImport(loadjson(fileToLoad, opt));

% scale - scaling of physical dimensions of region and actor coordinates.
%   Needed to mitigate patch display problems. Recommend that smallest
%   object (non-molecule) to plot has dimension of order 1, so a system
%   defined on the order of microns should have scale equal to 1e6.
scale = 1e6;

% customFigProp - structure of figure properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildFigureStruct for structure fields and their default values.
customFigProp = [];

% customAxesProp - structure of axes properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildAxesStruct for structure fields and their default values.
customAxesProp = [];

% regionToPlot - array of indices of regions to be plotted. No regions need
%   to be plotted.
regionToPlot = 1:config.numRegion;

% customRegionProp - structure of region properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildDispStruct for structure fields and their default values.
customRegionProp = [];

% actorToPlot - array of indices of actors to be plotted. Indexing matches
%   the actor list in the original config file and is independent of
%   whether an actor is active or passive. Actors listed here will have
%   their shapes plotted but NOT their molecules (the latter is indicated
%   by the argument "passiveActorToPlot". Actors are drawn after regions,
%   so if an actor is defined by region(s) then it will be drawn on top of
%   its region(s). No actors need to be plotted
actorToPlot = 1:config.numActor;

% customActorProp - structure of actor properties to change from AcCoRD
%   defaults. Can be passed as empty if no defaults are to be changed. See
%   accordBuildDispStruct for structure fields and their default values.
customActorProp = [];

% cameraAnchor - cell array of cell arrays defining anchor points for
%   the camera display. Can be passed as an empty cell array. Each anchor
%   point is a cell array defining a complete set of camera settings, in
%   the format {'CameraPosition', 'CameraTarget', 'CameraViewAngle',
%   'CameraUpVector'}. See MATLAB camera documentation for more details.
cameraAnchor = {};

%% Call accordPlotEnvironment
[hFig, hAxes] = accordEmptyEnvironment(fileToLoad, scale, ...
    customFigProp, customAxesProp, ...
    regionToPlot, customRegionProp, actorToPlot, customActorProp,...
    cameraAnchor);