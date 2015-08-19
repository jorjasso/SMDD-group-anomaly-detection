function k=RBF_kernel(X,Z,kernelParam)
% author: jorge
k=exp(-kernelParam*sqdistAll(X,Z));
