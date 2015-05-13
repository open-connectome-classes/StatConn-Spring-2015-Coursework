classdef SBM < handle
    %SBM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        adjMatrix;
        k;
        labels;
        bernBlock;
        bernEdge;
    end
    
    methods
        function this = SBM(a)
            this.adjMatrix = a;
        end
        
        function SetLabels(this, l)
            this.labels = l;
            this.k = max(l);
        end
        
        function SetK(this, k)
            this.k = k;
            this.labels = kmeans(this.adjMat, k);
        end
        
        
        
        function FitData(this)
            % Produce distributions among blocks(clusters)
            this.bernBlock = zeros(this.k);
            for i = 1 : this.k
                for j = 1 : this.k
                    subMat = this.adjMatrix(find(this.labels == i), find(this.labels == j));
                    this.bernBlock(i,j) = mean(subMat(:));
                end
            end
            
            % Produce distribution among edges
            this.bernEdge = zeros(size(this.adjMatrix));
            for i = 1 : this.k
                for j = 1 : this.k
                    this.bernEdge(find(this.labels == i), find(this.labels == j)) = this.bernBlock(i,j);
                end
            end
        end
        
        function ShowBernBlock(this)
            imagesc(this.bernBlock, [0, 1]);
            colormap gray
            axis equal tight
        end
        
        function ShowBernEdge(this)
            imagesc(this.bernEdge, [0, 1]);
            colormap gray
            axis equal tight
        end
        
        
        
        function adjMat = DrawAdjMat(this)
            adjMat = rand(size(this.bernEdge));
            for i = 1 : this.k
                for j = 1 : this.k
                    temp = adjMat(find(this.labels == i), find(this.labels == j));
                    temp(temp < this.bernBlock(i,j)) = 0;
                    temp(temp > 0) = 1;
                    adjMat(find(this.labels == i), find(this.labels == j)) = temp;
                end
            end
        end
        
        
        % The log-likelihood of an adjacency matrix(data) given the this SBM model
        function [ logLik, aic ] = LogLik(this, adjMat)
            edgeLik = adjMat .* this.bernEdge + (1-adjMat) .* (1-this.bernEdge);
            edgeLik(~edgeLik) = 0.000001;
            logLik = sum(log(edgeLik(:)));
            aic = -2 * logLik + 2 * this.k;
        end
        
        
        function vect = GetBlockFeature(this)
            indMat = triu(ones(size(this.bernBlock)));
            vect = this.bernBlock(find(indMat));
        end
    end
    
end

