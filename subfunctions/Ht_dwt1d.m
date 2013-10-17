% Written by Kun Qiu, ISU, April 2009

function alpha = Ht_dwt1d(At,b,h,L,m,n)
% alpha: sparse vector which is taken from Wavelet transform of image x
% A: random projection matrix K x N
% h: scaling filter
% L: level of decomposition
% m, n: size of image
% Return b: vector K x 1

% converting measurements into samples
if ~isa(At, 'function_handle')
    x=At*b;
else
    x=At(b);
end
% converting samples into wavelet coefficients (sparse representation)
im=reshape(x,m,n);
u=mdwt(im,h,L);
alpha=u(:);


