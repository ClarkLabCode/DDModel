# DDModel
Displacement Detector Model for LC11
## Instructions
Please clone or download the repository, and extend stimuli.zip to create "stimuli" folder directly under "DDModel" directory.
"model" folder has the function that runs numerical simulation of the displacement detector model (and some alternative models) reported in Tanaka & Clark (2020) Curr. Biol.
The scripts were created on Matlab 2019b.
"stimuli" folder should include folders for stimulus sets to be fed into the model. Each stimulus folder has a pair of xtPlot.xtp file and epochNames.mat file.
.xtp files are just large CSV with space-time plot of visual stimuli at 180Hz/5 degree resolution, flattened into 2D. Please read the comments in model/DDModel2DbyEpoch.m for more details. epochNames.mat contains a cell that has the names of the stimulus epochs included in the corresponding .xtp file.
