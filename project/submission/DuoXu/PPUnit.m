classdef PPUnit < handle
    %PROUNIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dFvF = [];
        
        pCorr = [];
        pCorrBin = [];
        
        pcVect = [];
        pcVal = [];
        pcThreshold = [];
        pcSignif = [];
        
        icSignal = [];
        icWeight = [];
    end
    
    methods
        function this = PPUnit(dfvf)
            if nargin > 0
                this.dFvF = dfvf;
            end
        end
        
        function ShowDfvf(this)
            dFluorPlotPopTraces(r, this.dFvF');
        end
        
        
        function PCorr(this)
        	this.pCorr = corr(this.dFvF);
            this.pCorrBin = zeros(size(this.pCorr));
            this.pCorrBin(this.pCorr > 0.5 | this.pCorr < -0.5) = 1;
        end
        
        function ShowCorr(this)
            imagesc(this.pCorr, [ -1 1 ]);
            axis equal tight
        end
        
        function ShowCorrBin(this)
            imagesc(this.pCorrBin, [ 0 1 ]);
            colormap gray
            axis equal tight
        end
        
        
        function PCA(this)
            % For concatenated trace
            [ this.pcVect, ~, this.pcVal ] = pca(this.dFvF);
            this.pcVal = this.pcVal';
            
            % Determine major PCs by Marcenko-Pastur distribution
            [ numRow, numCol ] = size(this.dFvF);
            q = numRow / numCol;
            this.pcThreshold = (1 + sqrt(1/q))^2;
            this.pcSignif = max(find(this.pcVal > this.pcThreshold));
        end
        
        function ShowEigVect(this, n)
            if nargin < 2
                n = length(this.pcVal);
            end
            imagesc(repmat(this.pcVal(1:n), length(this.pcVal), 1) .* this.pcVect(:,1:n));
        end
        
        function ShowEigVal(this, n)
            if nargin < 2
                n = length(this.pcVal);
            end
            plot(this.pcVal(1:n), 'o');
            hold on
            plot([ 0, n+1 ], [ this.pcThreshold, this.pcThreshold ]);
            hold off;
            xlim([ 0, n+1 ]);
        end
        
        
        function ICA(this, opt)
            if nargin < 2
                opt = 'none';
            end
            if strcmp(opt, 'reduce')
                [ this.icSignal, this.icWeight, ~ ] = ...
                    fastica(this.dFvF', 'lastEig', this.pcSignif, 'verbose', 'off');
            else
                [ this.icSignal, this.icWeight, ~ ] = ...
                    fastica(this.dFvF', 'numOfIC', this.pcSignif, 'verbose', 'off');
            end
        end
        
        function ShowIC(this)
            imagesc(this.icWeight);
        end
        
        
        function CleanRawData(this)
            this.dFvF = [];
%             this.pCorr = [];
        end
    end
    
end

