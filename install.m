%INSTALL Script to install logLIRA on first-time use.
installLogLIRA();

function installLogLIRA()
    fprintf('Installing logLIRA... ');

    requiredAddons = {'Uniform Manifold Approximation and Projection (UMAP)'};
    installedAddons = matlab.addons.installedAddons();
    for i = 1:numel(requiredAddons)
        if ~contains(requiredAddons{i}, installedAddons.Name)
            fprintf('Aborted.\n');
            error('logLIRA:install:missingAddOn', '%s not found. Install it via Matlab Add-On Manager.', requiredAddons{i});
        end
    end

    logLIRAPath = pwd();
    cd('..')

    if contains(logLIRAPath, strsplit(path, pathsep))
        fprintf('Skipped.\n');
        fprintf('logLIRA already installed.\n');
    else
        addpath(logLIRAPath);
        status = savepath();
        if status == 0
            fprintf('Done.\n');
        else
            fprintf('Aborted.\n');
            fprintf('Matlab path could not be saved.\n');
        end
    end

    cd(logLIRAPath);
end