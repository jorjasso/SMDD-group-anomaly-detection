function k=computeKernel(X,Z,kernelOp,kernelParam)
%            X,Z =  data matrices
%            kernelParam = kernel parammeter
%            kernelOp = 1 linear kernel, 2 polinomial kernel, 3 gaussian
%            kernel, 4 gaussian kernel with spherical normalization
% output    K= kernel gram matrix
% author: jorge.jorjasso@gmail.com
[L,Borrar]=size(X);
[LL,Borrar]=size(Z);

%muX=mean(X) or another estimator
%covX= cov(X) or another estimator

switch lower(kernelOp)
    case 1% linear kernel
        k=linear_kernel(X,Z);
    case 2 %polinomial kernel
        k=poly_kernel(X,Z,kernelParam);
    case {3,4} % gaussian kernel %option four for spherical normalization
        k=RBF_kernel(X,Z,kernelParam);
        
        %                       case 5 % Nystrom method in approximating the kernel  matrix (rbf kernel)
        %                          kernel = struct('type', 'rbf', 'para', 1/kernelParam);
        %                          k = INys(kernel,data, m, 'k');
        %                          k(1:10,1:10)
        %                          val=sum(sum(k))/(L*LL);
        %                          valTr=sum(sum(k))/(L*(LL-1));
        %                          trK=sum(diag(k))/(L-1)-valTr; %K_trace(i)=diag(kern)'*p'-G(i,i);
        
    otherwise
        disp('Unknown method.')
end
