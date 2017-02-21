# Support Measure Data Description:
A one-class classifier for group anomaly detection

## Synopsis
Experiments of the paper: Support Measure Data Description, A one-class classifier for group anomaly detection.

---
## Prerequisites
* SVM-KM from http://asi.insa-rouen.fr/enseignants/~arakoto/toolbox/  (included)
* MATLAB 2013a for the experiments
* R for the plots


## Code Example

Run the Matlab file main.m. That would generate <nowiki>*</nowiki>.mat within the /output folder.
Alternatively, you can use the screen program in linux as follows:



```matlab
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(1,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(2,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(3,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(4,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(5,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentSloan(1)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentSloan(2)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentSloan(3)"
```

## Generating a cvs file with the results
The following script will generate <nowiki>*</nowiki>.csv files within the /output folder with the results

```matlab
mat = dir('output/*.mat')
 
for q = 1:length(mat) 
    cont = load(strcat('output/',mat(q).name))
    dataName = mat(q).name
    strcat('output/',mat(q).name,'.csv')
    filename=strcat('output/',mat(q).name,'.csv')
    fid = fopen(filename, 'w') ;
    dlmwrite(filename,cont.statistics) 
    fclose(fid) ;
end
```
## Generating some  plots
The plots in the paper were generated with the R script
```R
plots.R
```

## Motivation

This project is part of my research on kernels on probability distributions

## Contributor
jorge.jorjasso@gmail.com

