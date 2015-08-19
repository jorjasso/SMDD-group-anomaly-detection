function k=poly_kernel(X,Z,kernelParam)
% author: jorge
k=((X*Z'+kernelParam(1)).^kernelParam(2));