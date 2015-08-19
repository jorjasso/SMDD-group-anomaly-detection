function [Krr,Kre,Kee,K_traceTraining,K_traceTest]=kernelPM(S, STest,kernelOp,kernelParam)
% kernel on probability measures
% % References: Support measure data description support measure description for group anomaly detection
% input  
%     S      = trainig set of replicates: S={S{1},S{2},...,S{N}}. Each S{i} has L replicates for observartion i
%     STest	     = test set of replicates
%     kernelOp = 1 lineal, 2 = polinomyal, 3 = RBF kernel, 4 = RBF kernel with spherical normalization
% output
% 	Krr = kernel matrix between training samples S
% 	Kre = kernel matrix between training samples and test samples
%	Kee = kernel matrix between test samples
%       K_traceTraining = trace of the covariance operator in the RKHS of training set
%	K_traceTest = trace of the covariance operator in the RKHS of test set
% jorge.jorjasso@gmail.com



%Krr=kernel for training samples
%----------------------
n=length(S);

K_trace=zeros(n,1);
Krr=zeros(n,n);
for i=1:n
    for j=i:n %the matrix is symetric
        [L,Borrar]=size(S{i}); [LL,Borrar]=size(S{j});
        k=computeKernel(S{i},S{j},kernelOp,kernelParam);
        Krr(i,j)= sum(sum(k))/(L*LL); 
        Krr(j,i)=Krr(i,j);
        if i==j       
             valTr=sum(sum(k))/(L*(LL-1));
             % trace of the covariance operator in the RKHS
             %_------------------------------------------
             trK=sum(diag(k))/(L-1)-valTr; %K_trace(i)=diag(kern)'*p'-G(i,i);
             K_traceTraining(i)=trK;
        end
    end
end

%Kre= kernel between training and test samples
%---------------------------------------------
m=length(STest);
Kre=zeros(n,m);
for i=1:n
    for j=1:m                
        [L,Borrar]=size(S{i}); [LL,Borrar]=size(STest{j});
        k=computeKernel(S{i},STest{j},kernelOp,kernelParam);
        Kre(i,j)= sum(sum(k))/(L*LL); 
    end
end


%Kee= kernel between test samples
%--------------------------------
k = cellfun(@(x) computeKernel(x,x,kernelOp,kernelParam),STest ,'UniformOutput', false );
L_vec = cell2mat(cellfun(@(x) size(x,1), STest ,'UniformOutput', false));
k_val=cell2mat(cellfun(@(x) sum(sum(x)), k ,'UniformOutput', false));
k_diag=cell2mat(cellfun(@(x) sum(diag(x)), k ,'UniformOutput', false));
Kee=k_val./(L_vec.^2);

% trace para el test set
%---------------------------------------
valTr=k_val./(L_vec.*(L_vec-1));
trK=k_diag./(L_vec-1)-valTr; %K_trace(i)=diag(kern)'*p'-G(i,i);
K_traceTest=   trK;        

%kernel with spherical normalization only for the RBF kernel
if (kernelOp==4)
    diag_Krr=diag(Krr);
    diag_Kee=Kee';
    Krr=Krr./sqrt(diag_Krr*diag_Krr');
    Kre=Kre./sqrt(diag_Krr*diag_Kee');
    Kee=Kee./sqrt(diag_Kee.*diag_Kee)';%only need the diagonal elements
end

