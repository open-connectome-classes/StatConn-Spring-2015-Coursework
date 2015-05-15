function [ corr] = ConVogHard_rQAP( A,B,m )

% [corr,iter] = ConVogHard( A,B,m ) is the syntax.
%  A,B are (m+n)x(m+n) adjacency matrices,
% loops/multiedges/directededges allowed.
% It is assumed that the first m vertices of A's graph
% correspond respectively to the first m vertices of B's graph,
% corr gives the vertex correspondences
% For example, corr=[ 1 2 7 16 30 ...
% means that the vtx1ofA-->vtx1ofB, 2-->2, 3-->7, 4-->16, 5-->30
%  example: EXECUTE the following:
% >> v=[ [1:5] 5+randperm(400)]; B=round(rand(405,405));A=B(v,v);
% >> [corr,P] = ConVogHard( A,B,5 ) ; [v; corr]
% ready June 1, 2012   (Donniell's code)


[totv,~]=size(A);
n=totv-m;

%The block submatrices of A and B
A12=A(1:m,m+1:m+n);
A21=A(m+1:m+n,1:m);
A22=A(m+1:m+n,m+1:m+n);
B12=B(1:m,m+1:m+n);
B21=B(m+1:m+n,1:m);
B22=B(m+1:m+n,m+1:m+n);

% Note
% AntiSeed Permutation Matrix = [0 x; y z]

%Maximum number of iterations
patience=1;

% change ratio in function value for stopping criterion
%If change ratio is lower, the iterations terminate
tol=1E-3;
epsilon=0.01;
P=ones(n,n)/n;
toggle=1;
iter=0;


while (toggle==1)&(iter<patience)
    iter=iter+1;
    %Compute Gradient at P
    Grad=A22*P*B22'+A22'*P*B22+A21*B21'+A12'*B12;
    
    %Compute Q (Q_tilde in paper) that maximizes trace(Q'Grad)
    ind=lapjv(-Grad,0.01);    
    Q=eye(n);
    Q=Q(ind,:);
    
    % Compute polynomial coefficients c, d ,e ,u and v
    c=trace(A22'*P*B22*P');
    d=trace(A22'*Q*B22*P')+trace(A22'*P*B22*Q');
    e=trace(A22'*Q*B22*Q');
    u=trace(P'*A21*B21'+P'*A12'*B12);
    v=trace(Q'*A21*B21'+Q'*A12'*B12);

    % Solve for alpha in the quadratic equation
    alpha=-(d-2*e+u-v)/(2*(c-d+e));
    
    %Function value at alpha=0
    f0=e+v;
    
    %current function value (value at alpha=1)
    f1=c+u;
    %Funcion value at alpha estimate
    falpha=(c-d+e)*alpha^2+(d-2*e+u-v)*alpha+e+v;
    
    
    newfval=0;
    % The next P estimate(P_tilde in paper) 
    if ((alpha>0)&&(falpha>f0)&&(falpha>f1))
        %The next P estimate is linear combination of P and Q
        P=alpha*P+(1-alpha)*Q;
        newfval=falpha;
    elseif (f0>(f1+epsilon))
        %The next P estimate is Q
        P=Q;
        newfval=f0;
    else
        %f1 is already a (local) maximum 
        % Terminate the FW algorithm
        toggle=0;
    end
    if ((abs(newfval-f1)/abs(f1))<tol)
        %change in function value is negligible
        % Terminate the FW algorithm
        toggle=0;
    end
end

%Project the doubly stochastic matrix to set of permutation matrices
%by solving linear assignment problem for the last time.
corr=lapjv(-P,0.01);

% The first m  vertices in A are already matched with the first m in B
% Add the matchings found by Frank Wolf
corr=[ 1:m,  m+corr];

