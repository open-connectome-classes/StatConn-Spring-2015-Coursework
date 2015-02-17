
clear all;clc;

a=[0 1 1 1 0 0 0 0 0; 1 0 1 0 1 0 1 1 1;1 1 0 0 0 1 0 0 0;1 0 0 0 0 0 1 0 0;0 1 0 0 0 0 0 1 0;0 0 1 0 0 0 0 0 1; 0 0 0 1 0 0 0 1 1; 0 0 0 0 1 0 1 0 1;0 0 0 0 0 1 1 1 0;];
temp=kmeans(a,3);

for i=1:5 %Say it failed if it outputs different categorization each time
    fail=(temp ~= kmeans(a,3));
    if sum(fail)~=9
        disp('fail')
    else
        disp('didn''t fail')
    end
end