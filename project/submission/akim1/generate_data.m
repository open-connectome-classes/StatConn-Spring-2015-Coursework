% akim1 150513
% statistical connectomics final project

n = 50; % size of grid
div = 4; % number of blocks (in a sbm)

p_mix = 0.1; % probability of adopting the random case over contiguous

p_aa = 0.1; % inter-block probability
p_ab = 0.001; % intra-block probability

% contiguous pixels
y_nonrand = zeros(n);
k = 1;
for i = 1:sqrt(div)
    for j = 1:sqrt(div)
            y_nonrand((i-1)*n/sqrt(div)+1:i*n/sqrt(div),...
                (j-1)*n/sqrt(div)+1:j*n/sqrt(div)) = k;
        k = k+1;
    end
end

% random pixels
y_rand = randi(div,n,n);
rand_mix = rand(length(y_rand));
y_mix = (rand_mix<p_mix).*y_rand+(rand_mix>=p_mix).*y_nonrand;
clear rand_mix;

% generate vector that contains block identity of each pixel
id_nonrand = reshape(y_nonrand, [numel(y_nonrand) 1]); % block assignments for non-random pixels
id_rand = reshape(y_rand, [numel(y_rand) 1]); % block assignments for random pixels
id_mix = reshape(y_mix, [numel(y_mix) 1]); % block assignments for mixture pixels


% generate a connectivity matrix
mat_nonrand = zeros(length(id_nonrand));
mat_rand = zeros(length(id_rand));
mat_mix = zeros(length(id_mix));


% populate connectivity matrix based on Bernoulli distribution
for i = 1:div
    
    % connectivity matrix for random pixels
    % inter-block
    temp_rand = triu(ones(size(mat_rand)),1).*...
        (rand(length(mat_rand))<p_aa);
    mat_rand(id_rand==i, id_rand==i) =...
        temp_rand(id_rand==i, id_rand==i);
    
    % intra-block
    temp_rand = triu(ones(size(mat_rand)),1).*...
        (rand(length(mat_rand))<p_ab);
    mat_rand(id_rand==i, id_rand~=i) =...
        temp_rand(id_rand==i, id_rand~=i);
    mat_rand = mat_rand | mat_rand';
    
    
    % connectivity matrix for contiguous pixels
    % inter-block
    temp_nonrand = triu(ones(size(mat_nonrand)),1).*...
        (rand(length(mat_nonrand))<p_aa);
    mat_nonrand(id_nonrand==i, id_nonrand==i) =...
        temp_nonrand(id_nonrand==i, id_nonrand==i);
    
    % intra-block
    temp_nonrand = triu(ones(size(mat_nonrand)),1).*...
        (rand(length(mat_nonrand))<p_ab);
    mat_nonrand(id_nonrand==i, id_nonrand~=i) =...
        temp_nonrand(id_nonrand==i, id_nonrand~=i);
    mat_nonrand = mat_nonrand | mat_nonrand';
    
   
    % connectivity matrix for contiguous pixels
    % inter-block
    temp_mix = triu(ones(size(mat_mix)),1).*...
        (rand(length(mat_mix))<p_aa);
    mat_mix(id_mix==i, id_mix==i) =...
        temp_mix(id_mix==i, id_mix==i);
    
    % intra-block
    temp_mix = triu(ones(size(mat_mix)),1).*...
        (rand(length(mat_mix))<p_ab);
    mat_mix(id_mix==i, id_mix~=i) =...
        temp_mix(id_mix==i, id_mix~=i);
    mat_mix = mat_mix | mat_mix';
end