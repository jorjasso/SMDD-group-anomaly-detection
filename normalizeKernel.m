function [Krr, Kre,Kee]=normalizeKernel(Krr, Kre, Kee)
% spherical normalization of data
% jorge.jorjasso@gmail.com
 diag_Krr=diag(Krr);
    diag_Kee=Kee';
    Krr=Krr./sqrt(diag_Krr*diag_Krr');
    Kre=Kre./sqrt(diag_Krr*diag_Kee');
    Kee=Kee./sqrt(diag_Kee.*diag_Kee)';%only need the diagonal elements

