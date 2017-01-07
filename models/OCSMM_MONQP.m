function [rho,alph1,Io,I1,I2,I3]=OCSMM_MONQP(G,S,C,kernelOp)
% % References: Support measure data description support measure description for group anomaly detection
% % References  One class support measure data description for group anomaly detection
%  (M4 in the kdd paper)
% input  
%     S      = set of replicates: S={S{1},S{2},...,S{N}}. Each S{i} has L replicates for observartion i
%     G	     = kernel gram matrix, between the set of replicates using a kernel on probability measures
%     C	     = regularization parameter
%     kernelOp = if kernelOp =1, then returns a value for c, and -1 otherwise
% output
% c=center,cNorm=norm of c, alph1 = alpha solutions of dual,
% Io={index for 0<alph<=C}, I1={index for
% 0<alph<C},I2={index for alph==C}, I3=ALLindex-I1UI2
% jorge.jorjasso@gmail.com
n=length(S);
%[~,G]=kernelPM(S,mu,Sigma,S,mu,Sigma,kernelOp,kernelParam);
% correction due to numerical error for the eigen values
G=G+eps*eye(n);

[~,p]=chol(G);
if (p~=0)
    for j=eps*2.^[0:1:40]
        G=G+j*eye(n);
        [~,p]=chol(G);
        if p==0
            break
        end
    end
end

e=ones(n,1);

%optimization
%-----------------
l = 10^-11;
verbose = 0;
 
D=zeros(n,1);
[alph, lambda, Io] = monqp(G,D,e,1,C,l,verbose);




%indices
%---------
th=getThreshold(alph);% getThreshold
%Io=find(alph>th);      % sv expansion 0<alph<=C
I2=find(alph>=C-th);  % sv errors    alph==C
I1=setdiff(Io,I2);     % sv  0<alph<C
I3=setdiff([1:1:n],Io);

if kernelOp==1  % kernel linear
    c=alph'*muTrain(Io,:);
else
    c=-1; %kernel no linear
end

%to standarize tue output
%------------------------
alph1=zeros(n,1);
alph1(Io)=alph;


rho=-lambda;
%pause



