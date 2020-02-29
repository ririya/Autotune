function addDirAndSubDirToPath(rootPath)

    addpath(rootPath);
    fileList = dir(rootPath);
    nFiles = length(fileList);
    
    for i=1:nFiles
       fileFullName = [rootPath '\' fileList(i).name];        
        if isdir(fileFullName)
            addpath(fileFullName);
            addDirAndSubDirToPath(fileFullName);
        end
        
    end

end

