%INSTALL Script to install logssar on first-time use.
fprintf('Installing logssar... ');

if contains(pwd(), strsplit(path, pathsep))
    fprintf('Skipped.\n');
    fprintf('logssar already installed.\n');
else
    addpath(pwd());
    fprintf('Done.\n');
end