function [] = loadData()
%Downloads data using wget

URL = 'https://www.dropbox.com/sh/idt3d0gylplyo31/AABrOJqzX07wQl7Ma3DZr6Zya/SandyaS72_data.zip?dl=0';
destination_folder = 'sandya_data';

stringCommand = ['wget ', '-O ', destination_folder, ' ', URL, ' --no-check-certificate'];
system(stringCommand);
unzip(destination_folder);

addpath('SandyaS72_data'); %<- This is the original folder name of what you zipped so it never changes

end

