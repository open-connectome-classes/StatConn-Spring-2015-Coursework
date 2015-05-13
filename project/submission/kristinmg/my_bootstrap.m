function pvalue = my_bootstrap(data1, data2,samplesize1,samplesize2,num_tests);


true_mean1 = mean(data1);
true_mean2 = mean(data2);
t_true = abs(true_mean1-true_mean2);

data = [data1; data2];
numsamples = samplesize1+samplesize2;
t = zeros(num_tests,size(data1,2));
for i = 1:num_tests
    % Randomly draw #numsamples number of samples (with replacement) from the data
    drawsamples = round((numsamples-1).*rand(numsamples,1)+1);

    samples_drawn = data([drawsamples],:);

    mean1 = mean(samples_drawn(1:samplesize1,:));

    mean2 = mean(samples_drawn(samplesize1+1:end,:));

    t(i,:) = abs(mean1-mean2);
end

pvalue = zeros(size(data1,2),1);
for j = 1:size(t,2)
    count = 0;
    for k = 1:size(t,1)
        if t(k,j) >= t_true(j)
            count = count+1;
        end
    end
    pvalue(j) = count/num_tests;
end

