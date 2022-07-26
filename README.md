# AFM-Line-Profile-Analyzer
All code should be run from the same folder: 
MFquestdlg.m is a 3rd party function that helps extend functionality towards UI features during accept/reject line profile calculations. 
AFMProfileAnalyzer is a helper function that is extended by BatchAFMProfileAnalyzer and Area_AFMProfile Analyzer. 
Batch(...) takes into account individual line profiles but accepts large amounts of files. 
Area_(...) does the same, but the input is 2N line profiles, where N is the number of individual nanosheets, and assumes that the data structure is comprised of pairs of line profiles for orthogonal lines over the same nanosheet.

Subsequent information on code usage and algorithm can be found in the pptx slidedeck. 
