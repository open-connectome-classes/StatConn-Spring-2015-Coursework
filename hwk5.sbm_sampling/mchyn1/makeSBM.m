function [ M , P ] = makeSBM( n,k )
%UNTITLED2 Summary of this function goes here
%   n = number verticies
%   k = number blocks

P = rand(k);
P = .5*(P+P'); %gives symmetric probability matrix
for i = 1:k
    for j = 1:k
        if i~=j
            P(i,j) = P(i,j)/3;
        elseif i==j && P(i,j)<.3
            P(i,i) = P(i,i)+.3;
        end
    end
end
v = rand(k,1);
for i = 1:k-1
    if v(i) < .2
        v(i) = v(i)+.2;
    elseif sum(v(1:k-1))>1
        v(:) = v(:)/2;
    end
end
v = round(v*n);
v(k) = n-sum(v(1:k-1));


m = cell(k,k);
county = 1;
for i = 1:k
    countx = 1;
    for j = 1:k
        m{i,j} = rand(v(county),v(countx))<P(i,j);
        countx = countx +1;
    end
    county = county+1;
end

M = cell2mat(m);



end

