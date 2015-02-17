%Michael Seung%
%Homework 1%

%generate adjacency matrix%
A = [1,1,0,0,1,0; 1,0,1,0,1,0; 0,1,0,1,0,0; 0,0,1,0,1,1; 1,1,0,1,0,0; 0,0,0,1,0,0];

%specify dimensions%
N = size(A,1);

%specifcy boundaries%
coords = [cos(2*pi*(1:N)/N); sin(2*pi*(1:N)/N)]';

%plot%
gplot(A, coords)
text(coords(:,1) - 0.1, coords(:,2) + 0.1, num2str((1:N)'), 'FontSize', 10)