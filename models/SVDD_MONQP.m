function [R,c,cNorm,alph1,Io,I1,I2,I3]=SVDD_MONQP(G,muTrain,C,kernelOp)
%        function [R,c,cNorm,alph,Io,I1,I2,I3]=SVDD_MONQP(G,muTrain,C,kernelOp)
%        solves SVDD problem with monqp solver
%input:  G       =  kernel matrix
%        muTrain = data matrix
%        C       = regularization parameter
%        kernelOp = kernel option (it is no necessary, only to return c if the linear kernel is used, other wise return -1)

[n,Borrar]=size(G);
e=ones(n,1);

G=G+0.000000000001*eye(n);

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

%optimization
%-----------------
l = 10^-11;
verbose = 0;
 
[alph, lambda, Io] = monqp(2*G,diag(G),e,1,C,l,verbose);


%representer theorem
%-----------
cNorm=alph'*G(Io,Io)*alph;
R=sqrt(lambda+cNorm);

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


