
function k=kernel(X,Z,kernelOp,kernelParam)
% return the kernel function
% input      X, Z = matrix of D-dimensional row samples 
%            kernelOp = 
%                       1 linear       
%                       2 polinomial   kernelParam = (d,gamma)  
%                       3 gaussian     kernelParam = (gamma)  
%            kernelOp = multivariate normal distribution
%                       4 linear       
%                       5 polinomial   kernelParam = (d,gamma)  
%                       6 gaussian     kernelParam = (gamma) 
% output    gram matrix G
[L,~]=size(X);
[LL,~]=size(Z);


switch lower(kernelOp)
    %-------------------------------------- 
          % p(x| mu,cov)=1/L, i.e uniform
          case 1% linear kernel
            %disp('linear kernel')
            k=X*Z';            
          case 2 %polinomial kernel
            k=((X*Z'+kernelParam(1)).^kernelParam(2));            
          case {3,4} % gaussian kernel %option four for spherical normalization
            %dX=diag(X*X'); dZ=diag(Z*Z'); 
            %XX=dX*ones(1,LL);ZZ=ones(L,1)*dZ';
            %k=exp(-kernelParam*(XX-2*X*Z'+ZZ));
             k=exp(-kernelParam*sqdistAll(X,Z)); 

            
          otherwise
            disp('Unknown method.')
        end
