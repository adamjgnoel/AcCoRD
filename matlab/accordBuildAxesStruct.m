function axesProp = accordBuildAxesStruct()
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordBuildAxesStruct.m - build structure with desired parameters for
%   plotting axes. All structure fields must be valid axes properties since
%   they will be directly applied to the axes handle via set. Field values
%   can be appended or changed within this function or by modifying the
%   output structure. The fields and parameters listed here are considered
%   "defaults".
%
% INPUTS
% NONE
%
% OUTPUTS
% axesProp - structure with common axes properties
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
axesProp.Box = 'on';        % If on, draw complete box outline along axis
                            % limits
axesProp.Clipping = 'on';   % If on, crop parts of curves that are defined
                            % beyond axis limits
axesProp.OuterPosition = [0 0 1 1]; % Relative position of axes within figure
axesProp.Visible = 'off';    % If off, plots will appear but the axis lines
                            % grid lines, labels, etc, will be invisible

% Font and text properties
axesProp.FontName = 'Helvetica';
axesProp.FontSize = 12;
axesProp.FontWeight = 'normal';
axesProp.TickLabelInterpreter = 'latex';

% Set range, scaling, and grid visibility along each axis
axesProp.XGrid = 'on';
axesProp.YGrid = 'on';
axesProp.ZGrid = 'on';
axesProp.XLim = [0 1];
axesProp.XLimMode = 'auto';     % Default is 'auto'. If 'manual', value of
                                % XLim will be used
axesProp.YLim = [0 1];
axesProp.YLimMode = 'auto';     % Default is 'auto'. If 'manual', value of
                                % YLim will be used
axesProp.ZLim = [0 1];
axesProp.ZLimMode = 'auto';     % Default is 'auto'. If 'manual', value of
                                % ZLim will be used
axesProp.XScale = 'linear';
axesProp.YScale = 'linear';
axesProp.ZScale = 'linear';

% Cameria view properties
axesProp.Projection = 'perspective';    % Default is 'orthographic'.
                                        % 'perspective' enables depth
                                        % perception
axesProp.CameraPosition = [0 0 0];
axesProp.CameraPositionMode = 'auto';   % Default is 'auto'
                                        % If 'manual', need to specify
                                        % CameraPosition
axesProp.CameraTarget = [0 1 0];
axesProp.CameraTargetMode = 'auto';     % Default is 'auto'
                                        % If 'manual', need to specify
                                        % CameraTarget
axesProp.CameraUpVector = [0 0 1];      % Default in 3D is [0 0 1]
axesProp.CameraUpVectorMode = 'auto';   % Default is 'auto'
                                        % If 'manual', need to specify
                                        % CameraUpVector
axesProp.CameraViewAngle = 6.6086;      % Default is 6.6086. Range [0, 180)
axesProp.CameraViewAngleMode = 'auto';  % Default is 'auto'
                                        % If 'manual', need to specify
                                        % CameraViewAngle