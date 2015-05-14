function [ P_test ] = Ptest( cell1,cell2, P_est_2toc1,P_est_2tou1,P_est_2to1,P_est_1to1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pairs=nchoosek(cell1+1,2);

P_2to1=P_est_2to1;
P_1to1=P_est_1to1;

y=0;
total=10000;

%Sampling a random graph===================================================
for tests=1:total
    connected=0;
    SBM_1block=rand(cell1)<P_1to1;
    SBM_12block=rand(cell1,cell2)<P_2to1;
    for c=1:cell1
        for r=1:c
            if r==c
                SBM_1block(r,c)=0;
            end
            SBM_1block(r,c)=SBM_1block(c,r);
            connected=connected+SBM_1block(r,c);
        end
    end

    %list off the pairs in terms of connected and unconnected
    cindex=0;
    uindex=connected;
    pairlist=zeros(pairs,2);
    for c=1:cell1
        for r=1:c        
            if(SBM_1block(r,c)==1)
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
            if(SBM_12block(pairlist(i,1),j))==1
                if(SBM_12block(pairlist(i,2),j))==1
                    pairSBM_12block(i,j)=1;
                end
            end
        end
    end

    in_cblock=pairSBM_12block(1:connected,:);
    in_ublock=pairSBM_12block((connected+1):pairs,:);

    P_2toc1=sum(sum(in_cblock))/(connected*cell2);
    P_2tou1=sum(sum(in_ublock))/((pairs-connected-1)*cell2);

    if (abs(P_2toc1-P_2tou1)>abs(P_est_2toc1-P_est_2tou1))
        y=y+1;
    end
    
end

P_test=y/total;

end

