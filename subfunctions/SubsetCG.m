function [y Residual i]=SubsetCG(x,y,A,Pt,IN,tol,verbose,MaxIt)
% Conjugate gradient algorithm to solve y=Dtmat * Dmat * x
%
% Dmat and Dtmat are either matrices, objects or function handles
% for Dmat matrices or object, Dtmat is ignored.
% type "help function_format" or "help object_format" for more information)
% y   - initial value
% IN  - indices of elements in y to be adapted 
% if  IN =[] or not specified, then IN = find(y)
% if  y=0 or not specified, y is zero vector
% tol - tolerance, norm(ChangeInResidual)<tol to stop
% verbos - 1 to plot progress (0 is default)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if          isa(A,'float')      P =@(z) A*z;  Pt =@(z) A'*z;
elseif      isobject(A)         P =@(z) A*z;  Pt =@(z) A'*z;
elseif      isa(A,'function_handle') 
    try
        if          isa(Pt,'function_handle'); P=A;
        else        error('If P is a function handle, Pt also needs to be a function handle. Exiting.'); end
    catch error('If P is a function handle, Pt needs to be specified. Exiting.'); end
else        error('P is of unsupported type. Use matrix, function_handle or object. Exiting.'); end


if isempty(y) && exist('IN','var')
    y=Pt(zeros(size(x)));
elseif isempty(y) && ~exist('IN','var')
    error('y can not be empty if IN is not specified!')
elseif ~isempty(y) && ~exist('IN','var')
   IN=find(y);
elseif ~isempty(y) && exist('IN','var')
    yy=y;
    y=Pt(zeros(size(x)));
    y(IN)=yy(IN);
end
N=length(y);


if nargin <6 || isempty(tol)
    tol=1e-6;
end

if nargin<7
    verbose = 0;
end

maxIter     = length(IN);
lx          = length(x);
sigsize     = x'*x/lx;
Residual    = x-P(y);
    
DR          = Pt(Residual);

pDDp        = zeros(maxIter,1);
pDDp(1)     = 1;
p           = zeros(N,1);
MASK        = zeros(N,1);
MASK(IN)    = 1;

if toc == 0 
    tic;
end
t=toc;
i=1;
done=0;
while ~done 
     
     if i==1
         p(IN)  = DR(IN);
         Dp     = P(p);
     else
         Pd     = P(DR.*MASK);
         b      = -Dp'*Pd/pDDp(i-1);
         p(IN)  = DR(IN) + b*p(IN);
         Dp     = Pd + b * Dp;
     end

     pDDp(i)    = Dp'*Dp;     
     a          = Residual'*Dp/(Dp'*Dp);
     y          = y+a*p;
     oR         = Residual;
     Residual   = Residual-a*Dp;
     DR         = Pt(Residual);
     
         if ((oR'*oR - Residual'*Residual)/(lx*sigsize) < tol) || i>=maxIter
                done=1;
         elseif verbose && toc-t>10
                display(sprintf('Conjugate gradient solver. Iteration: %i. --- stoping criterion: %1.3g ',i ,(oR'*oR - Residual'*Residual)/(lx*sigsize))) 
                t=toc;
         end
         if  exist('MaxIt','var')
             if i>=MaxIt
                 done=1;
             end
         end
     i=i+1;
end
i = i - 1;
%
% Change history:
% Iteration 1, set update direction to gradient.
% Bug fixed in conjugate gradient direction calculation
%

