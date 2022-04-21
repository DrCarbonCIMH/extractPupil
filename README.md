# extractPupil
Matlab functions to identify/extract the pupil parameters in a video file 

It uses a pulse-coupled neural network (PCNN) to identify and extract the pupil frame by frame
"Main functions" are:
- ms_preprocIMG: for a quick preprocessing of a frame
- ms_findPup_rb: identifying and parameter estimation of pupil
- ms_plotPup_rb: for plotting of found pupil

For an example run "pipeline_example.m" in the example folder

Keep in mind, that the PCNN is an iterative process which needs stopping criteria followed by chosing which iteration step was the "correct/best" one (see "bestCandidate" in ms_findPup_rb).
You probably have to adapt some parameters for your used hardware, video files, etc.

When using this code, please cite: 
[CIT]

The "pipeline function ExtractPupil_LW_4NatCom" used in the paper [CIT] can be found in "pipeline_usedInPaper"

Written at the Central Institute of Mental Health (CIMH), Mannheim, Germany
by Markus Sack, Robert Becker, and Laurens Winkelmeier; RG Translational Imaging, Department Neuroimaging