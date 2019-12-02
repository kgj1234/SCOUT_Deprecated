# SCOUT

SCOUT (Single-Cell SpatiOtemporal LongitUdinal Tracking) consists of two parts, spatial filters applied to a base optical recording extraction algorithm to improve ROI selection, and a cell tracking module in which we apply a predictor-corrector methodology to track neural cells extracted from optical recordings across multiple sessions.

The first module significantly decreases the false discovery rate when extracting neural activity from single recordings, without significantly reducing the number of detected neurons.

The cell tracking module performs well in regards to false discovery rate, and percent detected rate, when compared to similar methods such as cellReg.

These methods are described in "Robust Population Single Neuronal Calcium Signal Extraction Using SCOUT Allows for Longitudinal Analysis of Behavior-associated Neural Ensemble Dynamics”, which is currently under review by Nature Methods.

./Demos contains demos of both modules, on simulated recordings.
