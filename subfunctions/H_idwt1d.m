% Written by Kun Qiu, ISU, April 2009

function b=H_idwt1d(A,alpha,h,L,m,n)

% alpha: sparse vector which is taken from Wavelet transform of image x
% A: random projection matrix K x N
% h: scaling filter
% L: level of decomposition
% m, n: size of image
% Return b: vector K x 1

% converting wavelet coefficients (sparse representation) into samples
u=reshape(alpha,m,n);   
im=midwt(u,h,L);

% converting samples into measurements
x = im(:);    
if ~isa(A, 'function_handle')
    b = A*x;
else
    b = A(x);
end
    
   