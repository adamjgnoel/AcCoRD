function videoObj = accordInitializeVideo(videoName, videoFormat, propChange)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordInitializeVideo.m - initialize MATLAB video object to make a video
%   of AcCoRD simulation output.
%
% INPUTS
% videoName - video filename. File extension can be omitted and implied by
%   format
% videoFormat - string of video format. 'Motion JPEG AVI' is the MATLAB
%   default, but 'MPEG-4' is recommended for Windows since it creates much
%   smaller video files.
% propChange - structure of properties to modify from their default values
%   Any field must be a standard Matlab video property and be the correct
%   size. Fields that match any of the default values below will be applied
%   in the order below. Fields that do not match the defaults below will
%   be applied after the defaults and in the order defined. By defauly,
%   AcCoRD creates videos with 24 frames per second and at the highest
%   available quality.
%
% OUTPUTS
% videoObj - opened video object.
%
% Last revised for AcCoRD v0.6 (public beta, 2016-05-30)
%
% Revision history:
%
% Revision v0.6 (public beta, 2016-05-30)
% - Created file
%
% Created 2016-05-19

%% Create Video Object and Set Default Values

videoObj = VideoWriter(videoName, videoFormat);

videoObj.FrameRate = 24; % MATLAB Default 30
switch videoFormat
    case {'Archival', 'Uncompressed AVI', 'Indexed AVI', 'Grayscale AVI'}
       % Lossless formats cannot set quality value
    otherwise
        videoObj.Quality = 100; % MATLAB Default 75.
end

%% Make Specified Changes to Defaults
if ~isempty(propChange)
    propFields = fieldnames(propChange);
    numProp = numel(propFields);
    for i = 1:numProp
        if strcmp(propFields{i},'Quality')
            % Can only set quality to lossy formats
            switch videoFormat
                case {'Archival', 'Uncompressed AVI', 'Indexed AVI', 'Grayscale AVI'}
                   % Ignore quality setting
                otherwise
                    videoObj.(propFields{i}) = propChange.(propFields{i});
            end
        else
            videoObj.(propFields{i}) = propChange.(propFields{i});
        end
    end
end

%% Open Video
open(videoObj);