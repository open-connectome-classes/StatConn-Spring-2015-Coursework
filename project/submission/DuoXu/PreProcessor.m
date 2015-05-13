classdef PreProcessor < handle
    %PPTRIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rawTableFolder;
        rawTableName;
        rawTable;
        
        btypeInd = [ 9, 10, 11, 12 ];    % hit, miss, falseAlarm, correctRejection
        
        uSession;
        uBtypes;
        uTrials;
    end
    
    methods
        function this = PreProcessor()
            
        end
        
        
        function ImportData(this)
            if exist('dataPath.mat', 'file')
                load('dataPath.mat');
            else
                rtFolder = '';
            end
            
            [ this.rawTableName, this.rawTableFolder ] = uigetfile('*.mat', 'Import Imaging Data', rtFolder);
            load(fullfile(this.rawTableFolder, this.rawTableName));
            this.rawTable = D;
            
            rtFolder = this.rawTableFolder;
            save('dataPath.mat', 'rtFolder');
            
            FillingData(this);
        end
        
        
        function FillingData(this)
            % Initiate data structure
            this.uSession = PPUnit(this.rawTable.deltaFOverF);
            this.uBtypes = cell(length(this.btypeInd),1);
            this.uTrials = cell(length(this.btypeInd),1);
            
            % Looping through each behavioral type
            for i = 1 : length(this.btypeInd)
                % Extract sub-table of one behavioral type from rawTable
                ind = find(this.rawTable.(this.btypeInd(i)));
                currentBtype = this.rawTable(ind,:);
                numFrame = max(currentBtype.frameNumInTrial);
                numTrial = length(currentBtype.(4)) / numFrame;
                
                % Filling data into each data structure levels
                this.uBtypes{i} = PPUnit(currentBtype.deltaFOverF);
                this.uTrials{i} = cell(numTrial,1);
                for j = 1 : numTrial
                    this.uTrials{i}{j} = PPUnit(currentBtype.deltaFOverF((j-1)*numFrame+1 : j*numFrame, :));
                end
            end
        end
        
        
        function Calculate(this)
            this.uSession.PCorr();
%             this.uSession.PCA();
%             this.uSession.ICA();
            
            for i = 1 : length(this.uBtypes)
                this.uBtypes{i}.PCorr();
%                 this.uBtypes{i}.PCA();
%                 this.uBtypes{i}.ICA();
                
                for j = 1 : length(this.uTrials{i})
                    this.uTrials{i}{j}.PCorr();
%                     this.uTrials{i}{j}.PCA();
%                     this.uTrials{i}{j}.ICA();
                end
            end
        end
        
        
        function CleanRawData(this)
            this.uSession.CleanRawData();
            this.rawTable = [];
            
            for i = 1 : length(this.uBtypes)
                this.uBtypes{i}.CleanRawData();
                for j = 1 : length(this.uTrials{i})
                    this.uTrials{i}{j}.CleanRawData();
                end
            end
        end
        
        
        function ClearResult(this)
            this.uSession = [];
            this.uBtypes = [];
            this.uTrials = [];
        end
        
        
        function corrMats = GetBtypeCorrBins(this)
            corrMats = zeros([ size(this.uBtypes{1}.pCorr) length(this.uBtypes) ]);
            for i = 1 : size(corrMats,3)
                corrMats(:,:,i) = this.uBtypes{i}.pCorrBin;
            end
        end
        
    end
    
end

