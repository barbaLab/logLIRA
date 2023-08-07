%INSTALL Script to install logssar on first-time use.
fprintf('Installing logssar... ');

if contains(pwd(), strsplit(path, pathsep))
    fprintf('Skipped.\n');
    fprintf('logssar already installed.\n');
else
    addpath(pwd());
    status = savepath();
    if status == 0
        fprintf('Done.\n');
    else
        fprintf('Aborted.\n');
        fprintf('Matlab path could not be saved.\n');
    end
end