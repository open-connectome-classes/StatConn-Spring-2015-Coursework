

num_iterations = 5;
totSize = 400;
numSeeds = 0:num_iterations;


num_correct_seeds = [];
num_correct_antis = [];
for ns = numSeeds;
    
    cur_trial_numCorrect_seed = [];
    cur_trial_numCorrect_anti = [];
    
    for n = 1:5
        v= [1:totSize]; %[ [1:numSeeds] numSeeds+randperm(totSize-numSeeds)]; 
        B=round(rand(totSize,totSize));
        A=B(v,v);
        
        
        % seeds 
         [corr] = ConVogHard_rQAP( A,B,ns ); 
         cur_trial_numCorrect_seed = [cur_trial_numCorrect_seed sum(v == corr)];
        
        

        % antiseeds 
        cur_trial_numCorrect_anti = [];
        % change the permutation from an identity matrix to an flipped
        % diagonal to fix the antiseed block diagonal
        C = [B(:,ns+1:totSize) B(:,1:ns) ];
        anti_v = [[ns:totSize] [0:ns]];
        [corr] = AntiSeed( A,C,ns ); 
        
        cur_trial_numCorrect_anti = [cur_trial_numCorrect_anti sum(v  == corr)];

    end
    
    num_correct_seeds = [num_correct_seeds mean(cur_trial_numCorrect_seed)];
    num_correct_antis = [num_correct_antis mean(cur_trial_numCorrect_anti)];
end


figure
plot (numSeeds, num_correct_seeds/totSize, numSeeds, num_correct_antis/totSize)
title('Number of Seeds vs Accuracy')
xlabel('Number of Seeds/Anti-Seeds')
ylabel('Fraction of correctly matched Vertices')

legend('Seeded FAQ','Anti-Seeded FAQ','Location','northeast')

