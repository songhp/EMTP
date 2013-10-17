function [s_hat,A_index,Count,Golden_Iter,USS]=ADORE(H,Ht,y,varargin)
% Automatic Double Overrelaxation ECME (ADORE) Thresholding routine: 
% Coded by Kun Qiu (kqiu@iastate.edu)
% updated Feb 15, 2010
% 
% 
% Function usage
% =======================================================
% INPUT(compulsory):
% H:                              the sensing matrix (H can be either an matrix or function handle that computes H*x)
% Ht:                             the adjoint function handle that computes H^T*y (if H is a matrix, just input [])
% y:                               the measurement column vector
% 
% INPUT(optional):
% 'SearchLen':             the minimum length of the Golden search interval
%                                   (default=4)
% 'IsHOrthonormal':     whether H is orthonormal: 
%                                   0: no
%                                   1: yes
%                                   (default=1)
% 'InitialSig':                 the mx1 initial signal estimate
%                                   (default=zeros(m,1))
% 'Thresh':                   tolerance threshold for stopping the DORE algorithm
%                                   (default=1e-14)
% ========================================================
% OUTPUT:
% s_hat:              the signal estimate
% A_index:          the estimated support set
% Count:             Count of the number of EM iterations
% Golden_Iter:   the number of golden section iterations
% USS:               the value of the USS function
% ========================================================



if (nargin-length(varargin))~=3
    error('Missing required inputs!');
end

if isa(H, 'function_handle')
    N=length(y);
    m=length(Ht(y));
else
    [N,m]=size(H);
end


%Setting default values for the optional inputs
SearchLen=4;
IsHOrthonormal=1;
z_init=zeros(m,1);
thresh=1e-14;

%Read the optional inputs
if (rem(length(varargin),2)==1)
    error('Optional inputs must go by pairs!');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case upper('SearchLen')
                SearchLen=varargin{i+1};
                if SearchLen<4
                    SearchLen=4
                end
            case upper('IsHOrthonormal')
                IsHOrthonormal=varargin{i+1};
            case upper('InitialSig')
                z_init=varargin{i+1};
            case upper('Thresh')
                thresh=varargin{i+1};
            otherwise
                error(['Unrecognized optional input: ''' varargin{i} '''']);
        end
    end
end

ytyN=y'*y/N;

Count=0;
Golden_Iter=0;

Lower=0;
Upper=round(N/2);
alpha=(3-sqrt(5))/2;
Len=Upper-Lower;
mid1=round(Lower+alpha*Len);
mid2=round(Upper-alpha*Len);


USS_store(1)=-inf; s_hat_store{1}=[];  A_index_store{1}=[];
[s_hat_store{2},A_index_store{2},Count_temp,loglikelihood,delta2]=DORE(H,Ht,y,mid1,'IsHOrthonormal',IsHOrthonormal,'InitialSig',z_init,'Thresh',thresh,'visibility',0);
USS_store(2)=-0.5*(N-mid1-2)*log(delta2/ytyN)-0.5*mid1*log(N/m);
Golden_Iter=Golden_Iter+1;
Count=Count+Count_temp;
[s_hat_store{3},A_index_store{3},Count_temp,loglikelihood,delta2]=DORE(H,Ht,y,mid2,'IsHOrthonormal',IsHOrthonormal,'InitialSig',z_init,'Thresh',thresh,'visibility',0);
USS_store(3)=-0.5*(N-mid2-2)*log(delta2/ytyN)-0.5*mid2*log(N/m);
Golden_Iter=Golden_Iter+1;
Count=Count+Count_temp;
USS_store(4)=-inf; s_hat_store{4}=[];  A_index_store{4}=[];

while Len>SearchLen
    if USS_store(2)<USS_store(3)
        Lower=mid1; USS_store(1)=USS_store(2); s_hat_store{1}=s_hat_store{2};  A_index_store{1}=A_index_store{2};
        mid1=mid2;  USS_store(2)=USS_store(3);  s_hat_store{2}=s_hat_store{3};  A_index_store{2}=A_index_store{3};
        Len=Upper-Lower;
        mid2=round(Upper-alpha*Len);
        [s_hat_store{3},A_index_store{3},Count_temp,loglikelihood,delta2]=DORE(H,Ht,y,mid2,'IsHOrthonormal',IsHOrthonormal,'InitialSig',z_init,'Thresh',thresh,'visibility',0);
        USS_store(3)=-0.5*(N-mid2-2)*log(delta2/ytyN)-0.5*mid2*log(N/m);
        Golden_Iter=Golden_Iter+1;
        Count=Count+Count_temp;
    else
        Upper=mid2; USS_store(4)=USS_store(3); s_hat_store{4}=s_hat_store{3};  A_index_store{4}=A_index_store{3};
        mid2=mid1;  USS_store(3)=USS_store(2);  s_hat_store{3}=s_hat_store{2};  A_index_store{3}=A_index_store{2};
        Len=Upper-Lower;
        mid1=round(Lower+alpha*Len);
        [s_hat_store{2},A_index_store{2},Count_temp,loglikelihood,delta2]=DORE(H,Ht,y,mid1,'IsHOrthonormal',IsHOrthonormal,'InitialSig',z_init,'Thresh',thresh,'visibility',0);
        USS_store(2)=-0.5*(N-mid1-2)*log(delta2/ytyN)-0.5*mid1*log(N/m);
        Golden_Iter=Golden_Iter+1;
        Count=Count+Count_temp;
    end
end

[USS,ind]=max(USS_store);
s_hat=s_hat_store{ind};
A_index=A_index_store{ind};


