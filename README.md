# RHC_Echo_Optim - JPhysiol Manuscript Release
This code simulates patient number 27 from the manuscript. 
Simply run the SmithRed4_OptimScript.m file in Matlab and it will produce a plot of the optimized fit from the manuscript.

If the user wants to adjust parameters then the user must:
1) change the flags in the first column of the MP27_InputData.txt file according to the input file header text
2) adjust the parameter value in the SmithRed4_OptimScript.m file
3) rerun the script

If the user wants to input their own patient data then the user must:
1) input the general patient data in column 2 of the input text file
2) input the RHC data in column 3 of the input text file
3) input the TTE data in column 4 of the input text file
The user can then make the parameters adjustable as noted above
Also to note the delimiter between columns in the input file is a tab.

Please feel free to contact the authors if a full optimization code is needed.
