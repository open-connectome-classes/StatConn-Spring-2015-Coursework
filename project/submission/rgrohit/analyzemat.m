function [ P_est_2toc1,P_est_2tou1,P_est_2to1,P_est_1to1] = analyzemat( SBM,cell1,cell2 )
%Annalyzemat turns the SBM into a matrix that blocked based on whether
%cells are connected or not
%   Creates a matrix with a block of connected primary cells and a block of
%   unconnected primary cells. Using that, it finds estimated probability
%   of type 2 cell connecting to either set of type 1 pairs. Also finds
%   overall probability of type 2 cell connecting to pairs of type 1 pairs

%=======================================================================

%make a list of connected/unconnected pairs
connected=sum(sum(SBM(1:cell1,1:cell1)))/2;
pairs=nchoosek(cell1+1,2);
cindex=0;
uindex=connected;
pairlist=zeros(pairs,2);
for c=1:cell1
    for r=1:c        
        if(SBM(r,c)==1)
            cindex=cindex+1;
            pairlist(cindex,1)=r;
            pairlist(cindex,2)=c;
        else
            uindex=uindex+1;
            pairlist(uindex,1)=r;
            pairlist(uindex,2)=c;
        end
    end
end

%Make a pairSBM matrix that lets us differentiate the attatchment to connected PCs
%and unconnected PCs
pairSBM_12block=zeros(pairs,cell2);
for i=1:pairs
    for j=1:cell2
        if(SBM(pairlist(i,1),cell1+j))==1
            if(SBM(pairlist(i,2),cell1+j))==1
                pairSBM_12block(i,j)=1;
            end
        end
    end
end

%Get probability estimates
cpairs=(nchoosek(connected,2));
in_cblock=pairSBM_12block(1:connected,:);
in_ublock=pairSBM_12block((connected+1):pairs,:);

P_est_2toc1=sum(sum(in_cblock))/(connected*cell2); %estimate of probability of cell type 2 to connected pairs of cell type 1
P_est_2tou1=sum(sum(in_ublock))/((pairs-connected-1)*cell2); %estimate of probability of cell type 2 to unconnected pairs of cell type 1
P_est_2to1=sum(sum(SBM(1:cell1,(cell1+1):(cell1+cell2))))/(cell1*cell2); %estimate of probability of cell type 2 to cell type 1 overall
P_est_1to1=sum(sum(SBM(1:cell1,1:cell1)))/(cell1^2);

%figure of block of graph showing type 2 cell connection to pairs of type
%one
figure(2); hold all;
imagesc(pairSBM_12block);colormap('winter');
y=nchoosek(connected,2);x=linspace(1,cell2);plot(x,y,'k');
axis([1 cell2 1 pairs]);
ylabel('Pairs of excitatory cells'); xlabel('inhibitory cells');

end

