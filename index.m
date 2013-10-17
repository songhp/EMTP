%%
function index

% options: 'gm'=1 'gk'=2 'bm'=3 'bk'=4
% g: gaussian 
% b: binary 
% m: fixed measurements 
% k: fixed sparsity level
clear all
k = 30;
m = 80;
n = 400;
ensemble = 'USE';
opts.sigma = 0;
opts.gauss = 1;

% [x, y, A] = gen_signal(k, m, n, ensemble, opts);
[x, y, A] = gen_signal(k, m, n, ensemble, opts);

% 收敛误差
tol = 1e-6;

%  xt = GDE(A, y, k );

% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% z_init = zeros(n, 1);
% [xt,A_index,Count,GML_record] = AHT(H,Ht,y,z_init,1,1);
% Count
% SupportDetection(x, xt)

% tic
% Len_thresh=1;          
% z_init=zeros(m,1);
% [xt,A_index_ADORE,Count_ADORE,Golden_Iter,USS]=ADORE(A,[],y,'SearchLen',Len_thresh,'Thresh',tol,'IsHOrthonormal',0);
% SupportDetection(x, xt)
% toc
% 
% 
% tic
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% Len_thresh=1;          
% z_init=zeros(m,1);
% [xt,A_index_ADORE,Count_ADORE,Golden_Iter,USS]=ADORE(H,Ht,y,'SearchLen',Len_thresh,'Thresh',tol,'IsHOrthonormal',0);
% SupportDetection(x, xt)
% toc


% tic
% [xt, A_index,Count,loglikelihood,delta2] =  DORE(A, [], y, k, 'Thresh',tol,'visibility',0,'IsHOrthonormal',0);
% Count
% SupportDetection(x, xt)
% toc
% 
% tic
% invAAt      = inv(A*A');
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% invHHt      = @(z) invAAt*z;
% [xt, A_index,Count,loglikelihood,delta2] =  DORE(H, Ht, y, k, 'Thresh',tol,'Inverse',invHHt,'visibility',0,'IsHOrthonormal',0);
% Count
% SupportDetection(x, xt)
% toc
% 
% tic
% [xt, Out] = TEMP(A, y, 1, 'Tolerance', tol);
% Out.iter
% SupportDetection(x, xt)
% toc

% 
% tic
% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt, Out] = TEMP(H, y, 1, 'At', Ht, 'MaxIt',100);
% Out.iter
% SupportDetection(x, xt)
% toc

% tic
% [xt, Out] = EMTP(A, y, k, 1, 'Tolerance',tol);
% Out.iter
% SupportDetection(x, xt)
% toc

% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt, Out] = EMTP(H, y, k, 1, 'At', Ht, 'MaxIt',10);
% Out.iter
% SupportDetection(x, xt)


tic
[xt, Out] = EMTPbeta(A, y, 0.5, 'Tolerance',tol);
Out.iter
SupportDetection(x, xt)
toc

% H           = @(z) A*z;
% Ht          = @(z) A'*z;
% [xt, Out] = EMTPbeta(H, y, 0.5, 'At', Ht, 'MaxIt',100);
% Out.iter
% SupportDetection(x, xt)



% SupportDetection(x, xt)
% % SupportDetection(x, xbar) % 对比原始信号与重构信号



