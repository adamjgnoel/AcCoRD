function molStruct = accordBuildMolStruct(numActor, indToDisp, propChange)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordBuildMolStruct.m - build structure with parameters to define how
%   molecules should be plotted (as either markers or 3D objects)
%
% INPUTS
% numActor - number of actors whose observations are being recorded
% indToDisp - cell array of indices of molecule types that will be
%   displayed and thus need properties
% propChange - cell array of structure of properties to modify from their
%   default values. Any field must match one in molStruct and be the
%   same size.
%
% OUTPUTS
% molStruct - structure with marker display parameters
%
% Last revised for AcCoRD v1.2 (2018-05-30)
%
% Revision history:
%
% Revision v1.2 (2018-05-30)
% - renamed to accordBuildMolStruct from accordBuildMarkerStruct
% - added properties to enable drawing of molecules as 3D shapes as an
% alternative to drawing as plot markers. Added details for 'sphere'
%
% Revision v0.6 (public beta, 2016-05-30)
% - Created file
%
% Created 2016-05-18

molStruct = cell(1,numActor);

% Default box structure in case any molecules will be displayed as boxes
faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
vertices = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];

for i = 1:numActor
    %% Set Default Values

    indLength = length(indToDisp{i});

    molStruct{i} = struct('numToDisp', indLength, ...   % # of types being plotted
        'indToDisp', indToDisp{i}, ...                     % indices of those being plotted
        'dispStr', cell(1), ...                         % Display string (if legend turned on)
        'shape', cell(1), ...                           % Shape to draw molecules as {'marker', 'sphere'}
        'lineStyle', cell(1), ...                       % Style of edge lines for non-marker shapes
        'opaque', ones(1, indLength), ...               % Opaqueness of non-marker shapes
        'sphereSize', 15*ones(1, indLength), ...        % Number of points along dimension of sphere plot (only applies to sphere shape)
        'edgeColor', cell(1), ...                       % Molecule edge (border) color
        'faceColor', cell(1), ...                       % Molecule face color
        'marker', cell(1), ...                          % Marker shape (if shape is 'marker')
        'lineWidth', 0.5*ones(1, indLength), ...        % Molecule edge (border) width
        'size', 50*ones(1, indLength));                 % If shape is 'marker', marker size in points^2, where 1 point = 1/72 inch
                                                        % If shape is not 'marker', width in meters (scale will be applied)

    % Cell arrays need to be defined as a cells with cell arrays in order to
    % suppress the output structure being an array of structures
    molStruct{i}.dispStr{1} = cell(1,indLength);
    molStruct{i}.shape{1} = cell(1,indLength);
    molStruct{i}.lineStyle{1} = cell(1,indLength);
    molStruct{i}.edgeColor{1} = cell(1,indLength);
    molStruct{i}.faceColor{1} = cell(1,indLength);
    molStruct{i}.marker{1} = cell(1,indLength);

    for j = 1:indLength
        molStruct{i}.dispStr{1}{j} = sprintf('Molecule %d',indToDisp{i}(j));
        molStruct{i}.shape{1}{j} = 'marker';
        molStruct{i}.lineStyle{1}{j} = '-';
        molStruct{i}.edgeColor{1}{j} = 'black';
        molStruct{i}.faceColor{1}{j} = 'none';
        molStruct{i}.marker{1}{j} = 'o';
    end

    %% Make Specified Changes to Defaults
    if ~isempty(propChange{i})
        propFields = fieldnames(propChange{i});
        numProp = numel(propFields);
        for j = 1:numProp
            molStruct{i}.(propFields{j}) = propChange{i}.(propFields{j});
        end
    end
    
    %% Build shape coordinates for non-marker shapes
    for j = 1:indLength
        if strcmp(molStruct{i}.shape{1}{j}, 'sphere')
            [X, Y, Z] = sphere(molStruct{i}.sphereSize(j));
            X = molStruct{i}.size(j)/2*X;
            Y = molStruct{i}.size(j)/2*Y;
            Z = molStruct{i}.size(j)/2*Z;
            [curFaces, curVertices, ~] = surf2patch(X,Y,Z);
            molStruct{i}.faces{j} = curFaces;
            molStruct{i}.vertices{j} = curVertices;
            molStruct{i}.numVertices(j) = size(molStruct{i}.vertices{j},1);
        end
    end
end