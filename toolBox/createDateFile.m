function createDateFile(currentPath,fname)

% Create a new data folder if it doesn't exist
datafolder = [currentPath,'/',fname];
if exist(datafolder)==0
    mkdir(datafolder);
else
    disp('The data folder exists.');
end