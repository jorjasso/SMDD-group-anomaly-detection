function [R,c,cNorm,alph,Io,I1,I2,I3]=SMDDCCP(G,K_trace,S,kappa,C,kernelOp)
% References: Support measure data description support measure description for group anomaly detection
% support measure data description CCP approach (M3 in the kdd paper)
% input  
%     S      = set of replicates: S={S{1},S{2},...,S{N}}. Each S{i} has L replicates for observartion i
%    K_trace = vector (Nx1)of trace values of covariance matrices \Sigma_i
%    kappa   = vector of probabilistic bounds for the N constraints
%     G	     = kernel gram matrix, between the set of replicates using a kernel on probability measures
%     C	     = regularization parameter
%     kernelOp = if kernelOp =1, then returns a value for c, and -1 otherwise
% jorge.jorjasso@gmail.com
% output
% R= radius, c=center,cNorm=norm of c, alph, alpha solutions of dual,
% Io={index for 0<alph_ikappa_i<=C}, I1={index for
% 0<alph_ikappa_i<C},I2={index for alpha_ikappa_i=C}, I3=ALLindex-I1UI2
% 
n=length(S);
%[K_trace,G]=kernelPM(S,muTrain,Sigma,S,muTrain,Sigma,kernelOp,kernelParam);
% correction due to numerical error for the eigen values
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


%-----
l=chol(G); % Error =l'*l-K, gives alph'*l*l'*alph
e=ones(n,1);
%---optimization

cvx_begin
    %cvx_precision best
    cvx_quiet(true)
    variables  alph(n);  
    dual variable lambda
    minimize  (quad_over_lin(alph'*l',(alph'*e))-alph'*diag(G)-alph'* K_trace);
    subject to       
         zeros(n,1)<=alph.*kappa;
         kappa.*alph<=C*e;
    lambda: (alph.*kappa)'*e==1;     
cvx_end


th=getThreshold(alph);% getThreshold
Io=find(kappa.*alph>th);      % sv expansion 0<kappa.*alph<=C
I2=find(kappa.*alph>=C-th);          % sv errors    kappa.*alph==C
I1=setdiff(Io,I2);            % sv  0<kappa.*alph<C
I3=setdiff([1:1:n],Io);


%---Io, index set to compute the center
%Io=find((kappa.*alph)<=C&(kappa.*alph)>0.000009); % todos los vectores de soporte incluyendo alpha=C

% representer theorem for center c to comput ||c||^2
cNorm=(alph(Io))'*G(Io,Io)*(alph(Io))/sum(alph(Io))^2;
if kernelOp==1  % kernel linear %comprobar
    c=((alph(Io))'*muTrain(Io,:))/((alph(Io))'*e(Io));
else
    c=-1; %kernel no linear
end

%-----------------------
%RADIO in the RKHS
R=sqrt(-lambda);


%--------------------------------------------------------------------------
% %-----------------PRECISION NUMERICA
% for j=eps*2.^[40:-1:0]
%     I1=find((kappa.*alph)<(C-j)&(kappa.*alph)>j);
%     [s,t]=size(I1);
%     if s~=0
%         break
%     end
% end
% %----radius, I1, index set to compute the radius
% %var1=0;
% %compute all the radius for values in I1, the result is the mode
% %for i=1:length(I1)   
% %    var1(i) = sqrt((G(I1(i),I1(i))-2*alph(I1)'*G(I1,I1(i))/sum(alph(I1))+cNorm+K_trace(I1(i)))/kappa(I1(i)));
% %end
% %[a,in]=mode(round(var1*100));
% %R=var1(in)
% 
% %-------------I2, index to compute the errors
% I2=setdiff(Io,I1);  %is the same that I2=find(abs(alph-(C./kappa))<0.0009);
% % well estimated points
% I3=setdiff([1:1:n],Io);  % or [1:1:n]-I1 U I2
% 
% 
