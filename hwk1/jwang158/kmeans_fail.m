clear all;
close all;
clc;

RAND_SIZE = 10;
K = 4;

try
    for i = 1:1000
        A = randi([0 1],RAND_SIZE,RAND_SIZE);
        A = floor(0.5*(A + A')); % Symmetricalization
        c = kmeans(A,K);
    end
    fprinf('No A found that fails kmeans\n')
catch exception
    fprintf('A fails kmeans\n')
    A
end

x = linspace(0,2*pi,RAND_SIZE+1);
coords = [cos(x);sin(x)];
gplot(A,coords','o-black')
