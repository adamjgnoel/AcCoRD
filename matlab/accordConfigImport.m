function config =  accordConfigImport(configJSON)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordConfigImport.m - convert AcCoRD simulation configuration structure
%   that is imported via the JSON interface into a structure format that is
%   more consistent and easier to use. The converted structure is more
%   consistent in its use of cell arrays, has simpler field names, and adds
%   fields that may not have been explicitly included in the configuration
%   file (but which can be implied from the contents of the file)
%
% INPUTS
% configJSON - structure created by running "loadjson" function on a valid
%   AcCoRD configuration file
%
% OUTPUTS
% config - structure converted from configJSON. All content in configJSON
%   should be maintained. More fields are added whose values are implied by
%   the contents of configJSON (e.g., counts of number of regions and
%   actors)
%
% Last revised for AcCoRD v1.4 (2018-08-06)
%
% Revision history:
%
% Revision v1.4 (2018-08-06)
% - added A Priori Absorption as a surface reaction type. Also added
% corresponding unqiue parameters (surface reaction threshold type and
% threshold value)
% - fixed loading of surface reaction parameters when there are products
% that are not released away from the surface
%
% Revision v1.3 (2018-07-31)
% - fixed detection of actors that are defined as the union of more than 1
% region. JSON imported the list of regions as a string array and not a
% cell array
% - modified checks on some optional active actor parameters ("Slot
% Interval" and "Random Molecule Release Times?"), since they are only
% conditionally needed
%
% Revision v1.1 (2016-12-24)
% - added import of local (region) diffusion coefficients, and all flow
% parameters (global and local)
%
% Revision v1.0 (2016-10-31)
% - modified import of region placement parameters to accommodate simpler
% definitions allowed in config files
% - modified import of modulation parameters to accommodate addition of
% "burst" modulation, which does not modulate data
% - enable import of multiple diffusion coefficients for each molecule type
% (one per region)
% - added import of surface reaction diffusion coefficient when it has been
% specified
%
% Revision v0.7.0.1 (public beta, 2016-08-30)
% - corrected import of surface reaction probability type when the surface
% reaction type is normal
%
% Revision v0.6 (public beta, 2016-05-30)
% - Created file
%
% Created 2016-05-16

config.outputFilename = configJSON.Output_0x20_Filename;

config.numRegion = length(configJSON.Environment.Region_0x20_Specification);

%% Extract simulation control parameters
config.numRepeatPerSeed = configJSON.Simulation_0x20_Control.Number_0x20_of_0x20_Repeats;
config.dt = configJSON.Simulation_0x20_Control.Global_0x20_Microscopic_0x20_Time_0x20_Step;
config.tFinal = configJSON.Simulation_0x20_Control.Final_0x20_Simulation_0x20_Time;

%% Extract chemical properties
config.numMolTypes = ...
    configJSON.Chemical_0x20_Properties.Number_0x20_of_0x20_Molecule_0x20_Types;

config.diffusionCoeff = zeros(config.numRegion,config.numMolTypes);
for i = 1:config.numMolTypes
    config.diffusionCoeff(:,i) = ...
        configJSON.Chemical_0x20_Properties.Diffusion_0x20_Coefficients(i);
end

if isfield(configJSON.Chemical_0x20_Properties, 'Global_0x20_Flow_0x20_Vector')
    config.flowType = configJSON.Chemical_0x20_Properties.Global_0x20_Flow_0x20_Type;
else
    config.flowType = 'None';
end

config.flowVector = cell(config.numRegion,config.numMolTypes);
if ~strcmp(config.flowType, 'None') && isfield(configJSON.Chemical_0x20_Properties, 'Global_0x20_Flow_0x20_Vector')
    for i = 1:config.numRegion
        for j = 1:config.numMolTypes
            if isfield(configJSON.Chemical_0x20_Properties,'Does_0x20_Molecule_0x20_Type_0x20_Flow_0x3F_')
                if(configJSON.Chemical_0x20_Properties.Does_0x20_Molecule_0x20_Type_0x20_Flow_0x3F_(j))
                    config.flowVector{i,j} = configJSON.Chemical_0x20_Properties.Global_0x20_Flow_0x20_Vector;
                end
            else                
                config.flowVector{i,j} = configJSON.Chemical_0x20_Properties.Global_0x20_Flow_0x20_Vector;
            end
        end
    end
end

config.numChemRxn = ...
    length(configJSON.Chemical_0x20_Properties.Chemical_0x20_Reaction_0x20_Specification);
chemStruct = struct('label', '', 'bReversible', 0, 'revLabel', '', ...
    'bSurface', 0, 'surfRxnType', '', 'surfTransProb', '',...
    'bEverywhere', 0, 'numExceptions', 0, 'exceptionRegions', [], ...
    'reactants', zeros(1,config.numMolTypes), 'products', ...
    zeros(1,config.numMolTypes), 'bProdReleased', zeros(1,config.numMolTypes), ...
    'releaseType', '', 'rate', 0);
config.chemRxn = cell(1, config.numChemRxn);
for i = 1:config.numChemRxn
    config.chemRxn{i} = chemStruct;
    if iscell(configJSON.Chemical_0x20_Properties.Chemical_0x20_Reaction_0x20_Specification)
        curRxn = ...
            configJSON.Chemical_0x20_Properties.Chemical_0x20_Reaction_0x20_Specification{i};
    else
        curRxn = ...
            configJSON.Chemical_0x20_Properties.Chemical_0x20_Reaction_0x20_Specification(i);
    end
    
    config.chemRxn{i}.rate = curRxn.Reaction_0x20_Rate;
    
    config.chemRxn{i}.label = curRxn.Label;
    config.chemRxn{i}.bReversible = curRxn.Is_0x20_Reaction_0x20_Reversible_0x3F_;
    if config.chemRxn{i}.bReversible
        config.chemRxn{i}.revLabel = curRxn.Reverse_0x20_Reaction_0x20_Label;
    end
    
    config.chemRxn{i}.bSurface = curRxn.Surface_0x20_Reaction_0x3F_;
    if config.chemRxn{i}.bSurface
        config.chemRxn{i}.surfRxnType = ...
            curRxn.Surface_0x20_Reaction_0x20_Type;
        if ~strcmp(config.chemRxn{i}.surfRxnType, 'Normal')
            config.chemRxn{i}.surfTransProb = ...
                curRxn.Surface_0x20_Transition_0x20_Probability;
        end
        if isfield(curRxn, 'Surface_0x20_Reaction_0x20_Diffusion_0x20_Coefficient')
            config.chemRxn{i}.diffusion = ...
                curRxn.Surface_0x20_Reaction_0x20_Diffusion_0x20_Coefficient;
        end
    end
    
    config.chemRxn{i}.bEverywhere = curRxn.Default_0x20_Everywhere_0x3F_;
    if iscell(curRxn.Exception_0x20_Regions)
        config.chemRxn{i}.numExceptions = length(curRxn.Exception_0x20_Regions);
    elseif ~isempty(curRxn.Exception_0x20_Regions)
        config.chemRxn{i}.numExceptions = 1;
    else
        config.chemRxn{i}.numExceptions = 0;
    end
    
    config.chemRxn{i}.exceptionRegions = cell(1,config.chemRxn{i}.numExceptions);
    if config.chemRxn{i}.numExceptions == 1
        config.chemRxn{i}.exceptionRegions{1} = curRxn.Exception_0x20_Regions;
    else
        for j = 1:config.chemRxn{i}.numExceptions
            config.chemRxn{i}.exceptionRegions{j} = ...
                curRxn.Exception_0x20_Regions{j};
        end
    end
        
    for j = 1:config.numMolTypes
        config.chemRxn{i}.reactants(j) = curRxn.Reactants(j);
        config.chemRxn{i}.products(j) = curRxn.Products(j);
        if strcmp(config.chemRxn{i}.surfRxnType, 'Adsorbing') || ...
            strcmp(config.chemRxn{i}.surfRxnType, 'Desorbing') || ...
            strcmp(config.chemRxn{i}.surfRxnType, 'A Priori Absorbing')
            config.chemRxn{i}.bProdReleased(j) = ...
                curRxn.Products_0x20_Released_0x3F_(j);
        end
    end
    
    if (strcmp(config.chemRxn{i}.surfRxnType, 'Adsorbing') || ...
        strcmp(config.chemRxn{i}.surfRxnType, 'Desorbing') || ...
            strcmp(config.chemRxn{i}.surfRxnType, 'A Priori Absorbing')) && ...
        sum(config.chemRxn{i}.products) > 0 && ...
        sum(config.chemRxn{i}.bProdReleased(:)) > 0
        config.chemRxn{i}.releaseType =...
            curRxn.Release_0x20_Placement_0x20_Type;
    end
    
    if strcmp(config.chemRxn{i}.surfRxnType, 'A Priori Absorbing')
        config.chemRxn{i}.surfRxnThresholdType = ...
            curRxn.Surface_0x20_Reaction_0x20_Threshold_0x20_Type;
        config.chemRxn{i}.surfRxnThresholdValue = ...
            curRxn.Surface_0x20_Reaction_0x20_Threshold_0x20_Value;
    end
end

%% Extract Region Properties
config.subBaseSize = configJSON.Environment.Subvolume_0x20_Base_0x20_Size;

regionStruct = struct('label', '', 'parent', '', 'shape', '', ...
    'type', '', 'surfaceType', '', 'bMicro', 0,...
    'anchorCoor', zeros(1,3), 'subvolSizeInt', 0, ...
    'numSubDim', zeros(1,3), 'radius', 0);
config.region = cell(1, config.numRegion);
for i = 1:config.numRegion
    config.region{i} = regionStruct;
    if iscell(configJSON.Environment.Region_0x20_Specification)
        curRegion = configJSON.Environment.Region_0x20_Specification{i};
    else
        curRegion = configJSON.Environment.Region_0x20_Specification(i);
    end
    
    config.region{i}.label = curRegion.Label;
    config.region{i}.parent = curRegion.Parent_0x20_Label;
    config.region{i}.shape = curRegion.Shape;
    config.region{i}.type = curRegion.Type;
    if strcmp(config.region{i}.type, '3D Surface') || ...
        strcmp(config.region{i}.type, '2D Surface')
        config.region{i}.surfaceType = curRegion.Surface_0x20_Type;
    end
    
    if isfield(curRegion, 'Local_0x20_Diffusion_0x20_Coefficients')
        for j = 1:config.numMolTypes
            config.diffusionCoeff(i,j) = ...
                curRegion.Local_0x20_Diffusion_0x20_Coefficients(j);
        end
    end
    
    if isfield(curRegion, 'Local_0x20_Flow')
        numFlowExcept = length(curRegion.Local_0x20_Flow);
        for j = 1:numFlowExcept
            if numFlowExcept > 1
                curExcept = curRegion.Local_0x20_Flow{j};
            else
                curExcept = curRegion.Local_0x20_Flow;
            end
            for k = 1:config.numMolTypes
                if(curExcept.Is_0x20_Molecule_0x20_Type_0x20_Affected_0x3F_(k))
                    if strcmp(curExcept.Flow_0x20_Type, 'None')
                        config.flowVector{i,k} = [];
                    else
                        config.flowVector{i,k} = curExcept.Flow_0x20_Vector;
                    end
                end
            end
        end
    end
    
    if strcmp(config.region{i}.shape, 'Rectangle') || ...
        strcmp(config.region{i}.shape, 'Rectangular Box')
        config.region{i}.bMicro = curRegion.Is_0x20_Region_0x20_Microscopic_0x3F_;
        if isfield(curRegion, 'Number_0x20_of_0x20_Subvolumes_0x20_Per_0x20_Dimension')
            config.region{i}.numSubDim = zeros(1,3);
            for j = 1:3
                config.region{i}.numSubDim(j) = ...
                    curRegion.Number_0x20_of_0x20_Subvolumes_0x20_Per_0x20_Dimension(j);
            end
        else
            config.region{i}.numSubDim = ...
                [curRegion.Number_0x20_of_0x20_Subvolumes_0x20_Along_0x20_X ...
                curRegion.Number_0x20_of_0x20_Subvolumes_0x20_Along_0x20_Y ...
                curRegion.Number_0x20_of_0x20_Subvolumes_0x20_Along_0x20_Z];
        end
        config.region{i}.subvolSizeInt = curRegion.Integer_0x20_Subvolume_0x20_Size;
        config.region{i}.radius = 0;
    elseif strcmp(config.region{i}.shape, 'Sphere')
        config.region{i}.bMicro = true;
        config.region{i}.numSubDim = [1 1 1];
        config.region{i}.radius = curRegion.Radius;
    else
        warning('Region %d shape %s not recognized\n', i-1, config.region{i}.shape);
    end
    
    if isfield(curRegion, 'Anchor_0x20_Coordinate')
        config.region{i}.anchorCoor = zeros(1,3);
        for j = 1:3
            config.region{i}.anchorCoor(j) = ...
                curRegion.Anchor_0x20_Coordinate(j);
        end
    else
        config.region{i}.anchorCoor = [curRegion.Anchor_0x20_X_0x20_Coordinate ...
            curRegion.Anchor_0x20_Y_0x20_Coordinate ...
            curRegion.Anchor_0x20_Z_0x20_Coordinate];
    end
end

%% Extract Actor Properties
config.numActor = length(configJSON.Environment.Actor_0x20_Specification);
actorStruct = struct('bDefinedByRegions', 0, 'numRegion', 0, ...
    'regionList', [], 'shape', '', 'boundary', [], 'bActive', 0, ...
    'startTime', 0, 'bMaxActions', 0, 'maxActions', 0, 'bIndependent', 0, ...
    'actionInterval', 0, 'bRecord', 0, 'activeID', NaN, 'passiveID', NaN);
config.actor = cell(1, config.numActor);
config.numActive = 0;
config.numPassive = 0;
for i = 1:config.numActor
    config.actor{i} = actorStruct;
    if iscell(configJSON.Environment.Actor_0x20_Specification)
        curActor = configJSON.Environment.Actor_0x20_Specification{i};
    else
        curActor = configJSON.Environment.Actor_0x20_Specification(i);
    end
    
    % How is actor space defined
    config.actor{i}.bDefinedByRegions = ...
        curActor.Is_0x20_Location_0x20_Defined_0x20_by_0x20_Regions_0x3F_;
    if config.actor{i}.bDefinedByRegions
        % Actor is defined by a list of regions
        
        if isstring(curActor.List_0x20_of_0x20_Regions_0x20_Defining_0x20_Location)
            config.actor{i}.numRegion = ...
                length(curActor.List_0x20_of_0x20_Regions_0x20_Defining_0x20_Location);
        elseif ~isempty(curActor.List_0x20_of_0x20_Regions_0x20_Defining_0x20_Location)
            config.actor{i}.numRegion = 1;
        else
            config.actor{i}.numRegion = 0;
        end
        
        config.actor{i}.regionList = cell(1,config.actor{i}.numRegion);
        if config.actor{i}.numRegion == 1
            config.actor{i}.regionList{1} = ...
                curActor.List_0x20_of_0x20_Regions_0x20_Defining_0x20_Location;
        else
            for j = 1:config.actor{i}.numRegion
                config.actor{i}.regionList{j} = ...
                    curActor.List_0x20_of_0x20_Regions_0x20_Defining_0x20_Location{j};
            end
        end        
    else
        % Actor is defined by a shape
        config.actor{i}.shape = curActor.Shape;
        boundaryLength = length(curActor.Outer_0x20_Boundary);
        config.actor{i}.boundary = zeros(1, boundaryLength);
        for j = 1:boundaryLength
            config.actor{i}.boundary(j) = curActor.Outer_0x20_Boundary(j);
        end
    end
    
    config.actor{i}.bActive = curActor.Is_0x20_Actor_0x20_Active_0x3F_;
    if config.actor{i}.bActive
        config.numActive = config.numActive + 1;
        config.actor{i}.activeID = config.numActive;
    else
        config.numPassive = config.numPassive + 1;
        config.actor{i}.passiveID = config.numPassive;
    end
    
    config.actor{i}.startTime = curActor.Start_0x20_Time;
    config.actor{i}.bMaxActions = ...
        curActor.Is_0x20_There_0x20_Max_0x20_Number_0x20_of_0x20_Actions_0x3F_;
    if config.actor{i}.bMaxActions
        config.actor{i}.maxActions = ...
            curActor.Max_0x20_Number_0x20_of_0x20_Actions;
    end
    config.actor{i}.bIndependent = curActor.Is_0x20_Actor_0x20_Independent_0x3F_;
    config.actor{i}.actionInterval = curActor.Action_0x20_Interval;
    config.actor{i}.bRecord = ...
        curActor.Is_0x20_Actor_0x20_Activity_0x20_Recorded_0x3F_;
end

% Make second pass through actor list to store active and passive
% parameters
config.activeActor = cell(1, config.numActive);
config.passiveActor = cell(1, config.numPassive);
activeStruct = struct('bRandNumMolecules', 0, 'bRandReleaseTimes', 0, ...
    'releaseInterval', 0, 'slotInterval', 0, 'bBitsRandom', 0, 'probBitOne', 0, ...
    'bitSequence', [], 'numBits', 0, ...
    'modScheme', '', 'numModBits', 1, 'modStrength', 0, 'bReleaseType', ...
    zeros(1,config.numMolTypes), 'actorID', 0);
passiveStruct = struct('bRecordTime', 0, 'bRecordMolCount', zeros(1,config.numMolTypes), ...
    'bRecordMolPosition', zeros(1, config.numMolTypes), 'actorID', 0);
curActive = 1;
curPassive = 1;
for i = 1:config.numActor
    if iscell(configJSON.Environment.Actor_0x20_Specification)
        curActor = configJSON.Environment.Actor_0x20_Specification{i};
    else
        curActor = configJSON.Environment.Actor_0x20_Specification(i);
    end
    
    if config.actor{i}.bActive
        config.activeActor{curActive} = activeStruct;
        config.activeActor{curActive}.actorID = i;
        
        config.activeActor{curActive}.bRandNumMolecules = ...
            curActor.Random_0x20_Number_0x20_of_0x20_Molecules_0x3F_;
        if config.activeActor{curActive}.bRandNumMolecules
            config.activeActor{curActive}.bRandReleaseTimes = ...
                curActor.Random_0x20_Molecule_0x20_Release_0x20_Times_0x3F_;
        else
            config.activeActor{curActive}.bRandNumMolecules = false;
        end
        config.activeActor{curActive}.releaseInterval = ...
            curActor.Release_0x20_Interval;
        if config.activeActor{curActive}.bRandNumMolecules
            config.activeActor{curActive}.slotInterval = 0;
        else
            config.activeActor{curActive}.slotInterval = ...
                curActor.Slot_0x20_Interval;
        end
        
        config.activeActor{curActive}.modScheme = ...
            curActor.Modulation_0x20_Scheme;
        switch config.activeActor{curActive}.modScheme
            case 'Burst'
                % "Burst" modulation does not modulate bits
            otherwise
                config.activeActor{curActive}.bBitsRandom = ...
                    curActor.Bits_0x20_Random_0x3F_;
                if config.activeActor{curActive}.bBitsRandom
                    config.activeActor{curActive}.probBitOne = ...
                        curActor.Probability_0x20_of_0x20_Bit_0x20_1;
                else
                    config.activeActor{curActive}.numBits = ...
                        length(curActor.Bit_0x20_Sequence);
                    config.activeActor{curActive}.bitSequence = ...
                        zeros(1, config.activeActor{curActive}.numBits);
                    for j = 1:config.activeActor{curActive}.numBits
                        config.activeActor{curActive}.bitSequence(j) = ...
                            curActor.Bit_0x20_Sequence(j);
                    end
                end
                config.activeActor{curActive}.numModBits = ...
                    curActor.Modulation_0x20_Bits;
        end
        
        
        config.activeActor{curActive}.modStrength = ...
            curActor.Modulation_0x20_Strength;
        
        for j = 1:config.numMolTypes
            config.activeActor{curActive}.bReleaseType(j) = ...
                curActor.Is_0x20_Molecule_0x20_Type_0x20_Released_0x3F_(j);
        end
        
        curActive = curActive + 1;
    else
        config.passiveActor{curPassive} = passiveStruct;
        config.passiveActor{curPassive}.actorID = i;
        
        config.passiveActor{curPassive}.bRecordTime = ...
            curActor.Is_0x20_Time_0x20_Recorded_0x20_with_0x20_Activity_0x3F_;
        
        for j = 1:config.numMolTypes
            config.passiveActor{curPassive}.bRecordMolCount(j) = ...
                curActor.Is_0x20_Molecule_0x20_Type_0x20_Observed_0x3F_(j);            
            config.passiveActor{curPassive}.bRecordMolPosition(j) = ...
                curActor.Is_0x20_Molecule_0x20_Position_0x20_Observed_0x3F_(j);
        end
        
        curPassive = curPassive + 1;
    end
end