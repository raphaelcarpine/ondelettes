function update_successful = test()
%TEST Summary of this function goes here
%   Detailed explanation goes here

folderPath = '_divers/test_update';

update_successful = true;
try
    rmdir folderPath s
catch
    update_successful = false;
end

disp('success!');


end

