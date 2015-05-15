function [ corr ] = AntiSeed( A,B,m )
%UNTITLED Summary of this function goes here
%   Elan Hourticolon-Retzler 5/11/15
%   Based on code by Donniell

% m - number of anti seeds 

[totv,~]=size(A);
n=totv-m;

%The block submatrices of A and B
A11=A(1:m,1:m);
A12=A(1:m,m+1:m+n);
A21=A(m+1:m+n,1:m);
A22=A(m+1:m+n,m+1:m+n);

B11=B(1:m,1:m);
B12=B(1:m,m+1:m+n);
B21=B(m+1:m+n,1:m);
B22=B(m+1:m+n,m+1:m+n);

% Note
% AntiSeed Permutation Matrix = [0 x; y z]

% T = A*P*B'*P'
% trace(T) = tr(A12*y*B21'*x' + (A11*x+A12*z)*B22'*x') + tr(A22*y*(B11'*y'+b21'*z')+(A21*x+A22*z)*(B12'*y'+B22'*z') )
%
%          = tr( A12*y*B21'*x' ) 
%          + tr( A11*x*B22'*x' ) 
%          + tr( A12*z*B22'*x' ) 
%          + tr( A22*y*B11'*y' )
%          + tr( A22*y*b21'*z' )
%          + tr( A21*x*B12'*y' )
%          + tr( A21*x*B22'*z' )
%          + tr( A22*z*B12'*y' )
%          + tr( A22*z*B22'*z' )

%dT_dx = d/dx(  tr( A12*y*B21'*x' ) 
%             + tr( A11*x*B22'*x' ) 
%             + tr( A12*z*B22'*x' ) 
%             + tr( A21*x*B12'*y' ) 
%             + tr( A21*x*B22'*z' ) 
%            ) 
%
%    grad_x = (A12*y*B21') + (A11'*x*B22 + A11*x*B22') + (A12*z*B22') +
%    (A21'*y*B12)+(A21'*z*B22)



%dT_dy = d/dy(  tr( A12*y*B21'*x' )  ==> tr(B21*x*A12'*y')
%             + tr( A22*y*B11'*y' )
%             + tr( A22*y*B21'*z' )  ==> tr(B21*z*A22'*y')
%             + tr( A21*x*B12'*y' )
%             + tr( A22*z*B12'*y' ) 
%            )
%
% grad_y = (A12'*x*B21) + (A22'*y*B11 + A22*y*B11') + (A22'*z*B21) + (A21*x*B12' + A22*z*B12');


%dT_dz = d/dz(  tr( A12*z*B22'*x' )  ==> tr( B22*x*A12'*z' )
%             + tr( A22*y*B21'*z' )
%             + tr( A21*x*B22'*z' )
%             + tr( A22*z*B12'*y' )  ==> tr( B12*y*B22'*z' )
%             + tr( A22*z*B22'*z' )
%            )
%
%grad_z = (A12'*x*B22) + (A22*y*B21') + (A21*x*B22') + (A22*y*B12) + (A22'*z*B22 + A22*z*B22');

%Maximum number of iterations
patience=1;

% change ratio in function value for stopping criterion
%If change ratio is lower, the iterations terminate
tol=1E-3;
epsilon=0.01;
x=ones(m,n)/totv;
y=ones(n,m)/totv;
z=ones(n,n)/totv;
P=[zeros(m,m) x; y z];
toggle=1;
iter=0;


while (toggle==1)&(iter<patience)
    iter=iter+1;
    
    % extract x,y,z
    x = P(1:m,m+1:m+n);
    y = P(m+1:m+n,1:m);
    
    %Compute Gradient at P
%     q = A12*y*B21';
%     size(q)
%     w = (A11'*x*B22 + A11*x*B22');
%     size(w)
%     e = A12*z*B22';
%     size(e)
%     r = B12*y*A21';
%     size(r)
%     t = A21'*y*B12 + A21'*z*B22;
%     size(t)
% %     y =
   
    grad_x = (A12*y*B21') + (A11'*x*B22 + A11*x*B22') + (A12*z*B22') + (A21'*y*B12)+(A21'*z*B22);
    grad_y = (A12'*x*B21) + (A22'*y*B11 + A22*y*B11') + (A22'*z*B21) + (A21*x*B12' + A22*z*B12');
    grad_z = (A12'*x*B22) + (A22*y*B21') + (A21*x*B22') + (A22*y*B12) + (A22'*z*B22 + A22*z*B22');
   
%     grad_x = A12*y*B21'+()+A12*z*B22'+(A21'*y*B12)+(A21'*z*B22)
    
%     size(grad_x)
%     size(grad_y)
%     size(grad_z)
    Grad= [[zeros(m,m);y] [grad_x; grad_z]];
    
    %Compute Q (Q_tilde in paper) that maximizes trace(Q'Grad)
    ind=lapjv(-Grad,0.01);    
    Q=eye(totv); %changed from n to totv
    Q=Q(ind,:);
    
    % Compute polynomial coefficients c, d ,e ,u and v
%     c=trace(A22'*P*B22*P');
%     d=trace(A22'*Q*B22*P')+trace(A22'*P*B22*Q');
%     e=trace(A22'*Q*B22*Q');
%     u=trace(P'*A21*B21'+P'*A12'*B12);
%     v=trace(Q'*A21*B21'+Q'*A12'*B12);

    c=trace(A'*P*B*P');
    d=trace(A'*Q*B*P')+trace(A'*P*B*Q');
    e=trace(A'*Q*B*Q');
    u=trace(P'*A*B'+P'*A'*B);
    v=trace(Q'*A*B'+Q'*A'*B);

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
