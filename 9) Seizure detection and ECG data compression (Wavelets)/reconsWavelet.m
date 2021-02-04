function [recons,coeff_total] = reconsWavelet(vectorC,vectorL,M)

%%INPUTS
%M = Number of highest coefficients to be computed
%vectorC = Output of wavelet decomposition command, wavelet decomposition vector c
%vectorL = Output of wavelet decomposition command, bookkeeping vector l

%%OUTPUT
%recons = reconstructed signal 
%coeff_total = vector of all wavelet coefficients (approx and detail)

%--------------------------------------------------
A5 = appcoef(vectorC,vectorL,'db4',5);
[D5,D4,D3,D2,D1] = detcoef(vectorC,vectorL,[5,4,3,2,1]);

coeff_total = [A5;D5;D4;D3;D2;D1];

[~,coeff_ind] = sort(coeff_total,'descend','ComparisonMethod','abs');

coeff_total(coeff_ind(M+1:end))=0; 

recons = waverec(coeff_total,vectorL,'db4');

end