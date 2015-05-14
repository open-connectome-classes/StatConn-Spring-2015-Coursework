% akim1 150513

function ret = coarsen_conn(n_x, n_y, conn_mat)

new_matrix = zeros(n_x*n_y/4);

l = 1;
for j = 1:2:n_x
   
    for i = 1:2:n_y
        % figure out the mask in the pixel domain
        mask = zeros(n_y, n_x);
        mask(i:i+1,j:j+1) = ones(2, 2);
        mask_v = reshape(mask, [numel(mask) 1]);
        
        % determine the connectivity
        a = conn_mat(mask_v==1,:);
        
        % traverse and populate the new coarsened connectivity matrix
        inds = [find(a(1,:)) find(a(2,:)) find(a(3,:)) find(a(4,:))];
        inds = unique(inds);
        for k = inds
            % remap the connectivity to the original grid
            mask = zeros(n_y*n_x, 1);
            mask(k) = 1;
            mask = reshape(mask, n_y, n_x);
            
            % find the coordinate in the original grid
            [r, c] = find(mask);
            
            % calculate the index of the connectivity matrix of the
            % coarsened grid
            new_k = n_x/2*(ceil(c/2)-1)+ceil(r/2);
            
            % assign the connection
            new_matrix(l, new_k) = 1;
            new_matrix(new_k, l) = 1;
        end
        l = l+1;
    end
end
ret = new_matrix;