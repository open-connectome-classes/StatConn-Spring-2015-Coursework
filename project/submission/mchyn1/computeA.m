A=cell(40,1);
for i = 1:40;
    er=0;
    Y=M{i}'; %M is the the with 0s discarded
    Y=normc(Y);
    
    [tp, v] = size(Y);
    Spec = vgxset('n',v,'nAR',1);
    %nAR is 1 because we're only doing t-1, so Y starts at 2 and Y0 is the first point
    try
        [EstSpec, EstStdErrors, LLF, W] = vgxvarx(Spec,Y(2:tp,:),[],Y(1,:));
    catch
        er=1;
        disp('too many parameters for data, but still continuing');
    end
    %vgxdisp(EstSpec, EstStdErrors)
    
    %Since Y is a single path, W is the inferred path of n-dim.
    EstW = vgxinfer(EstSpec, Y(2:tp,:), [],Y(1,:));
    
    %EstSpec.AR{1} is the autoregression matrix that we can use as A :)
    if er == 0;
        A{i} = EstSpec.AR{1};
    end
end
