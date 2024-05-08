# `logLIRA`: LOGarithmic Linear Interpolation for the Rejection of Artifacts

`logLIRA` is a novel stimulus artifacts rejection algorithm designed specifically for electrophysiological recordings. This algorithm employs a piece-wise linear interpolation with logarithmically distributed interpolation points to effectively estimate the stimulus artifacts templates and decouple the stimulus artifacts and the underlying neural activity.

Thanks to its embedded mitigation of secondary artifacts step, `logLIRA` is especially useful in the recording and analysis of short-latency evoked activity.

## Requirements

To successfully install and run `logLIRA`, you need:
-   MATLAB (R2017a or higher)
-   [Uniform Manifold Approximation and Projection (UMAP)](https://it.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap/) add-on installed in MATLAB

## Installation

To use logLIRA in your electrophysiological data analysis pipeline, simply clone this repository and incorporate the algorithm into your codebase via the following steps:

1. Clone this git repository with the following command:
   ```
   git clone https://github.com/barbaLab/logLIRA.git
   ```
2. Open MATLAB and run the `install.m` script. If prompted, allow MATLAB to update its path.

## Usage

```matlab
% Load electrophysiological data and sampling frequency
load('.../data.mat', 'data', 'sampleRate');

% Load the electrical stimuli onsets
load('.../stimuli.mat', 'stimIdxs');

% Apply logLIRA artifact rejection
output = logLIRA(data, stimIdxs, sampleRate);

% Proceed with your analysis steps
...
```

> For more details about `logLIRA` usage, read the documentation via the `help` command directly from MATLAB or check out the `logLIRA.m` file.

## Citation

## Contacts

## License