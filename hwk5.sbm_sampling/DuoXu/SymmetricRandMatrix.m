function m = SymmetricRandMatrix(n)
%SYMMETRICRANDMATRIX Summary of this function goes here
%   Detailed explanation goes here

mDir = triu(rand(n));
m = mDir + tril(mDir',-1);


end

