%% EXPERIMENTS GMM DATA SET 
% Requiered: CVX toolbox from http://cvxr.com/cvx/
run /Users/jorgeluisguevaradiaz/Documents/GITProjects/cvx/cvx_startup.m
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
% Requiered:  MONQP  from http://asi.insa-rouen.fr/enseignants/~arakoto/toolbox/
addpath ./SVM-KM/
addpath ./datasets/
addpath ./experiments/
addpath ./models/
addpath ./kernels/
addpath ./utils/




%Distribution-based group anomaly detection
%------------
%(perform nRuns (200 in the paper) to get statistics;
nRuns=2;
experimentGroupAnomalyDetectionGMM(0,nRuns) 

%Point-based group anomaly detection
%--------------
experimentGroupAnomalyDetectionGMM(1,nRuns)


%Distribution-based group anomaly detection (Another)
%------------

experimentGroupAnomalyDetectionGMM(1,nRuns) 

%% EXPERIMENTS SLOAN DATA SET

experimentSloan(0,nRuns)
experimentSloan(1,nRuns)
