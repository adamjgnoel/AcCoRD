function [data, config] = accord_import(fileName, seedRange, bWrite)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accord_import.m - load simulation output from the AcCoRD simulator
% 					into MATLAB
%
% INPUTS
% fileName - "root" filename of simulation output, without the "_SEED##" suffix
% 			or "summary_" prefix. File location is relative
% 			to the AcCoRD "results" directory
% seedRange - array of integers listing the seed values to read in. Generally, the
% 			independent realizations of a single simulation can be stored in multiple
% 			files, and each file has a "_SEED##" suffix. The results from all files
% 			specified are aggregated here into a single set of output parameters.
%
% OUTPUTS
% data - structure with simulation parameters and output. Outputs are also saved to the
% 		mat-file [data.configName]_out.mat, where data.configName is printed up until the
%		first "."
% config - structure with simulation configuration defined in user configuration file
%
% Last revised for AcCoRD LATEST_VERSION
%
% Revision history:
%
% Revision LATEST_RELEASE
% - added bWrite input argument to control whether output is saved to
% mat-file
%
% Revision v0.4.1
% - changed output filename convention and location for increased
% flexibility
%
% Revision v0.3.1
% - header added
%
% Created 2015-05-09

% Add subdirectory with JSONlab files to path (needed for loadjson)
addpath('JSONlab');

data.numSeed = length(seedRange);
data.seed = seedRange;

%% Import information from summary file
% Summary file is in JSON format
opt.SimplifyCell = 1;
opt.FastArrayParser = 1;
opt.ShowProgress = 0;

summary = loadjson([fileName ...
    '_SEED' num2str(seedRange(1)) '_summary.txt'], opt);
data.configName = summary{1}.ConfigFile; % Original configuration file

% "Search" for configuration file
bFoundConfig = false;
if exist(data.configName, 'file')
    configFile = data.configName;
	bFoundConfig = true;
elseif exist(['config/' data.configName], 'file')
    configFile = ['config/' data.configName];
	bFoundConfig = true;
elseif exist(['../config/' data.configName], 'file')
    configFile = ['../config/' data.configName];
	bFoundConfig = true;
elseif exist(['../' data.configName], 'file')
    configFile = ['../' data.configName];
	bFoundConfig = true;
elseif exist(['../../' data.configName], 'file')
    configFile = ['../../' data.configName];
	bFoundConfig = true;
elseif exist(['../../config' data.configName], 'file')
    configFile = ['../../config' data.configName];
	bFoundConfig = true;
else
	warning(sprintf('The original configuration file %s used to run the simulation cannot be found', data.configName));
end

if bFoundConfig
	config = loadjson(configFile, opt);
else
	config = [];
end
%data.seed = summary{1}.SEED;
data.startTime = summary{1}.StartTime;
data.endTime = summary{2}.EndTime;

data.numRepeatSingle = summary{1}.NumRepeat;
data.numRepeat = summary{1}.NumRepeat * data.numSeed;

data.numActive = summary{2}.NumberActiveActor;
data.activeID = zeros(1,data.numActive);
data.activeMaxBits = zeros(1,data.numActive); % Max # of bits in any realization
activeBitStr = cell(1,data.numActive);
data.activeBits = cell(1,data.numActive);
for i = 1:data.numActive
    data.activeID(i) = summary{2}.ActiveInfo(i).ID;
    data.activeMaxBits(i) = summary{2}.ActiveInfo(i).MaxBitLength;
    activeBitStr{i} = repmat('%d',1,data.activeMaxBits(i));
    data.activeBits{i} = zeros(data.numRepeat, data.activeMaxBits(i));
end

data.numPassiveRecord = summary{2}.NumberPassiveRecord;
data.passiveRecordID = zeros(1,data.numPassiveRecord);
data.passiveRecordBTime = zeros(1,data.numPassiveRecord);
data.passiveRecordMaxCountLength = zeros(1,data.numPassiveRecord);
passiveDoubleStr = cell(1,data.numPassiveRecord);
passiveCountStr = cell(1,data.numPassiveRecord);
data.passiveRecordNumMolType = zeros(1,data.numPassiveRecord);
data.passiveRecordMolID = cell(1,data.numPassiveRecord);
data.passiveRecordBPos = cell(1,data.numPassiveRecord);
data.passiveRecordTime = cell(1,data.numPassiveRecord);
data.passiveRecordCount = cell(1,data.numPassiveRecord);
data.passiveRecordPos = cell(1,data.numPassiveRecord);
for i = 1:data.numPassiveRecord
    data.passiveRecordID(i) = summary{2}.RecordInfo(i).ID;
    data.passiveRecordBTime(i) = summary{2}.RecordInfo(i).bRecordTime;
    data.passiveRecordMaxCountLength(i) = summary{2}.RecordInfo(i).MaxCountLength;
    passiveDoubleStr{i} = repmat('%f',1,data.passiveRecordMaxCountLength(i));
    passiveCountStr{i} = repmat('%u64',1,data.passiveRecordMaxCountLength(i));
    data.passiveRecordNumMolType(i) = summary{2}.RecordInfo(i).NumMolTypeObs;
    data.passiveRecordMolID{i} = zeros(1, data.passiveRecordNumMolType(i));
    data.passiveRecordBPos{i} = zeros(1, data.passiveRecordNumMolType(i));
    if data.passiveRecordBTime(i)
        data.passiveRecordTime{i} = zeros(data.numRepeat,...
            data.passiveRecordMaxCountLength(i));
    end
    data.passiveRecordCount{i} = zeros(data.numRepeat,data.passiveRecordNumMolType(i),...
        data.passiveRecordMaxCountLength(i));
    data.passiveRecordPos{i} = cell(1,data.passiveRecordNumMolType(i));
    for j = 1:data.passiveRecordNumMolType(i)
        data.passiveRecordMolID{i}(j) = summary{2}.RecordInfo(i).MolObsID(j);
        data.passiveRecordBPos{i}(j) = summary{2}.RecordInfo(i).bRecordPos(j);
        if data.passiveRecordBPos{i}(j)
            data.passiveRecordPos{i}{j} = cell(data.numRepeat,data.passiveRecordMaxCountLength(i));
        end
    end
end

%% Initialize Matlab data structure based on summary data

%% Read in realization data
curReal = 1;
for s = 1:data.numSeed
    seed = seedRange(s);
    fid = fopen([fileName '_SEED' num2str(seed) '.txt']);
    for r = 1:data.numRepeatSingle
        % Read in realization line
        textscan(fid, '%*[^\n]', 1);
        for i = 1:data.numActive
            % Read in active actor label line
            textscan(fid, '%*[^\n]', 1);
            % Read in actor bits
            content = textscan(fid, activeBitStr{i}, 1, 'CollectOutput',1);
            data.activeBits{i}(curReal,:) = content{:};
        end

        for i = 1:data.numPassiveRecord
            % Scan in passive actor label line
            textscan(fid, '%*[^\n]', 1);
            if data.passiveRecordBTime(i)
                % Time is being recorded. Read in time label
                textscan(fid, '%*[^\n]', 1);
                % Read in timestamps
                content = textscan(fid, passiveDoubleStr{i}, 1, 'CollectOutput',1);
                data.passiveRecordTime{i}(curReal,:) = content{:};
            end
            for j = 1:data.passiveRecordNumMolType(i)
                % Read in molecule type line AND count label line
                textscan(fid, '%*[^\n]', 2);
                % Read in molecule counts
                content = textscan(fid, passiveCountStr{i}, 1, 'CollectOutput',1);
                data.passiveRecordCount{i}(curReal,j,:) = content{:};
                if data.passiveRecordBPos{i}(j)
                    % Positions are being recorded. Read in Position label
                    textscan(fid, '%*[^\n]', 1);
                    for k = 1:length(data.passiveRecordCount{i}(curReal,j,:))
                        data.passiveRecordPos{i}{j}{curReal,k} = zeros(data.passiveRecordCount{i}(curReal,j,k),3);
                        % Read in opening round bracket associated with
                        % current count
                        textscan(fid, '%*[(]', 1);
                        for l = 1:data.passiveRecordCount{i}(curReal,j,k)
                            % Scan to start of coordinate
                            textscan(fid, '%*[(]', 1);
                            % Read in coordinates of current molecule
                            content = textscan(fid, '%f', 3, 'CollectOutput',1,'Delimiter',',');
                            data.passiveRecordPos{i}{j}{curReal,k}(l,:) = content{:};
                            % Scan to start of coordinate
                            textscan(fid, '%*[)]', 1);
                        end
                        % Scan in next newline
                        textscan(fid, '%*[^\n]', 1);
                    end
                end
            end
        end
        curReal = curReal + 1;
    end
    fclose(fid);
end

%% Write output to a .mat file
if bWrite
    [~,configName,~] = fileparts(data.configName);
    save([configName '_out'], 'data', 'config');
end