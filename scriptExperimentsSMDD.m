%% EXPERIMENTS GMM DATA SET 
% Requiered: CVX toolbox from http://cvxr.com/cvx/
run   /home/jorjasso/cvx/cvx_startup.m
%run /home/jorjasso/Downloads/cvx/cvx_startup.m
% Requiered:  MONQP  from http://asi.insa-rouen.fr/enseignants/~arakoto/toolbox/
addpath ./SVM-KM/


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

experimentGroupAnomalyDetectionGMM(2,nRuns) 

%% EXPERIMENTS SLOAN DATA SET

experimentSloan(0,nRuns)
experimentSloan(1,nRuns)
