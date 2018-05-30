function [data, config] = accordImport(fileName, seedRange, bWrite)
%
% The AcCoRD Simulator
% (Actor-based Communication via Reaction-Diffusion)
%
% Copyright 2016 Adam Noel. All rights reserved.
% 
% For license details, read LICENSE.txt in the root AcCoRD directory
% For user documentation, read README.txt in the root AcCoRD directory
%
% accordImport.m - load simulation output from the AcCoRD simulator
% 					into MATLAB
%
% INPUTS
% fileName - "root" filename of simulation output, without the "_SEED##" suffix
% 			or "_SEED##_summary" suffix. File path must be included if the
% 			output is not in the current directory.
% seedRange - array of integers listing the seed values to read in. Generally, the
% 			independent realizations of a single simulation can be stored in multiple
% 			files, and each file has a "_SEED##" suffix. The results from all files
% 			specified are aggregated here into a single set of output parameters.
% bWrite - boolean to generate an output file or not.
%
% OUTPUTS
% data - structure with simulation parameters and output. Outputs are also saved to the
% 		mat-file fileName_out.mat, where fileName is printed up until the
%		first "."
% config - structure with simulation configuration defined in user configuration file
%
% Last revised for AcCoRD v1.2 (2018-05-30)
%
% Revision history:
%
% Revision v1.2 (2018-05-30)
% - updated code that writes to output file so that filename is not cropped
% if "fileName" has a period in it
%
% Revision v1.1 (2016-12-24)
% - changed import algorithm to read entire output file in one call and
% store contents in a cell array, and then scan the elements of the cell
% array. Now about an order of magnitude faster at reading the simulation
% output.
%
% Revision v1.0 (2016-10-31)
% - corrected opening comments about where the output file(s) need to be
%
% Revision v0.6 (public beta, 2016-05-30)
% - renamed from accord_import for consistency
% - added conversion of default config structure, as loaded from the JSON
% interpreter, into a stucture that is easier to work with
% - changed format of output filename to match output file instead of
% configuration
% - made output of active actor data sequence a user option. Accounted for
% by changing how data.numActive is written
% - removed call to sprintf within call to warning
%
% Revision v0.5 (2016-04-15)
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
	warning('The original configuration file %s used to run the simulation cannot be found', data.configName);
end

if bFoundConfig
	config = accordConfigImport(loadjson(configFile, opt));
else
	config = [];
end
%data.seed = summary{1}.SEED;
data.startTime = summary{1}.StartTime;
data.endTime = summary{2}.EndTime;

data.numRepeatSingle = summary{1}.NumRepeat;
data.numRepeat = summary{1}.NumRepeat * data.numSeed;

data.numActive = summary{2}.NumberActiveRecord;
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
    fileText = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    curRow = 1;
    for r = 1:data.numRepeatSingle
        % Read in realization line
        curRow = curRow + 1;
        for i = 1:data.numActive
            % Skip active actor label line
            curRow = curRow + 1;
            % Read in actor bits
            content = textscan(fileText{1}{curRow}, activeBitStr{i}, 1, 'CollectOutput',1);
            data.activeBits{i}(curReal,:) = content{:}';
            curRow = curRow + 1;
        end

        for i = 1:data.numPassiveRecord
            % Scan in passive actor label line
            curRow = curRow + 1;
            if data.passiveRecordBTime(i)
                % Time is being recorded. Read in time label
                curRow = curRow + 1;
                % Read in timestamps
                content = textscan(fileText{1}{curRow}, passiveDoubleStr{i}, 1, 'CollectOutput',1);
                data.passiveRecordTime{i}(curReal,:) = content{:};
                curRow = curRow + 1;
            end
            for j = 1:data.passiveRecordNumMolType(i)
                % Skip molecule type line AND count label line
                curRow = curRow + 2;
                % Read in molecule counts
                content = textscan(fileText{1}{curRow}, passiveCountStr{i}, 1, 'CollectOutput',1);
                data.passiveRecordCount{i}(curReal,j,:) = content{:};
                curRow = curRow + 1;
                if data.passiveRecordBPos{i}(j)
                    % Positions are being recorded. Skip Position label
                    curRow = curRow + 1;
                    for k = 1:length(data.passiveRecordCount{i}(curReal,j,:))
                        curPos = 1;
                        data.passiveRecordPos{i}{j}{curReal,k} = zeros(data.passiveRecordCount{i}(curReal,j,k),3);
                        for l = 1:data.passiveRecordCount{i}(curReal,j,k)
                            % Go to start of coordinate
                            curPos = curPos + 2;
                            % Read in coordinates of current molecule
                            [content,pos] = textscan(fileText{1}{curRow}(curPos:end), '%f', 3, 'CollectOutput',1,'Delimiter',',');
                            curPos = curPos + pos;
                            data.passiveRecordPos{i}{j}{curReal,k}(l,:) = content{:}';
                            % Move to end of coordinate
                            curPos = curPos + 1;
                        end
                        % Move to next line
                        curRow = curRow + 1;
                    end
                end
            end
        end
        curReal = curReal + 1;
        curRow = curRow + 1;
    end
end

%% Write output to a .mat file
if bWrite
    [~,fileName,fileEnd] = fileparts(fileName);
    if ~isempty(fileEnd)
        fileName = [fileName fileEnd];
    end
    save([fileName '_out.mat'], 'data', 'config');
end