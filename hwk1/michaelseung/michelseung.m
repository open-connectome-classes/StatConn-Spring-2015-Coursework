%Michael Seung%
%Homework 1%

A = [1,1,0,0,1,0; 1,0,1,0,1,0; 0,1,0,1,0,0; 0,0,1,0,1,1; 1,1,0,1,0,0; 0,0,0,1,0,0];
%Adjancy matrix%

N = size(A,1);
%define dimensions%

coords = [cos(2*pi*(1:N)/N); sin(2*pi*(1:N)/N)]';
%define bound%

gplot(A, coords)

text(coords(:,1) - 0.1, coords(:,2) + 0.1, num2str((1:N)'), 'FontSize', 10)
