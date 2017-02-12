%% EXPERIMENTS GMM DATA SET 
% Requiered: CVX toolbox from http://cvxr.com/cvx/
% Requiered:  MONQP  from http://asi.insa-rouen.fr/enseignants/~arakoto/toolbox/

run /Users/jorgeluisguevaradiaz/Documents/GITProjects/cvx/cvx_startup.m
addpath ./SVM-KM/
addpath ./datasets/
addpath ./experiments/
addpath ./models/
addpath ./kernels/
addpath ./utils/
addpath ./output/

%Point-based group anomaly detection
%--------------
N=300;
percentAnomalies=10;
experimentsGADGMM(1,percentAnomalies,N)
experimentsGADGMM(2,percentAnomalies,N)
experimentsGADGMM(3,percentAnomalies,N)
experimentsGADGMM(4,percentAnomalies,N)


%Distribution-based group anomaly detection
%------------
N=300;
percentAnomalies=10;
experimentsGADGMM(5,percentAnomalies,N)


%% EXPERIMENTS SLOAN DATA SET

experimentSloan(1)
experimentSloan(2)
experimentSloan(3)

%% EXPERIMENTS USING SCREEN in linux
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(1,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(2,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(3,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(4,10,300)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentsGADGMM(5,10,300)"

screen  -d -m matlab -nodisplay -nosplash -r  "experimentSloan(1)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentSloan(2)"
screen  -d -m matlab -nodisplay -nosplash -r  "experimentSloan(3)"

%% write results of experiments into a csv files

for q = 1:length(mat) 
    cont = load(strcat('output/',mat(q).name))
    dataName = mat(q).name
    strcat('output/',mat(q).name,'.csv')
    filename=strcat('output/',mat(q).name,'.csv')
    fid = fopen(filename, 'w') ;
    dlmwrite(filename,cont.statistics) 
    fclose(fid) ;
end

%% Remark:
% the values of kappa can ne changed in the file getModelSetup.m (if options)


%%

%%%%%
% max frequency processors: cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq
Experiments GMM 
                                    11:30           14:00       18:26         09:18         10:08
-tolstoi(3800000),   (2,1) 3    1 | (2,2) 4  1  | (2,6) 2  1 | (2,10) 3  1 | (3,10) 1  1 | (5,2)  2  1 | (8,8)  1  1 
-hulk1(2301000),     (1,1) 5    2 | (1,2) 5  2  | (1,6) 2  2 | (1,9)  4  2 | (2,7)  4  2 | (3,11) 5  2 | (7,1)  4  2 
-bool (3900000),     (1,1) 2    4 | (1,1) 5  4  | (1,3) 2  4 | (1,6) 2   4 | (1,12) 5  4 | (2,9)  2  4 | (3,16) 3  2 
-escience4 (4400000) (1,1) 3    3 | (1,2) 3  3  | (1,4) 5  3 | (1,9) 2   3 | (2,6)  1  3 | (3,10) 5  3 | 
-puchkin                                                                   
     
Experiments Sloan          9:25         10:23
-dostoievski           (3,12) 1    1 |(5,16) 5    1   (finished)
-puchkin               (1,13) 5    2 |(2,9) 5     2 | (3,3) 5  3
