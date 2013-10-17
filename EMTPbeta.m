function  [x,Out] = EMTPbeta(A, y, beta, varargin)
% iterative detection-estimation
% 

%%%%%%%%%% Input specification %%%%%%%%%%%%%%%%%%%
%  A:	measurement matrix or operator denoted by a structure.
%  y:	measured data (could be noisy)
%  opts (optional):
%	opts.x0:	initial signal
%   opts.xs:    true signal; if provided, debug flag (Dflag) can be set to true
%	opts.Th:	threshold
%   opts.maxit: maximal number of iteration
%
%%%%%%%%%% Input specification %%%%%%%%%%%%%%%%
%
%  x:	reconstructed signal
%  Out: (optional)
%     Out.iter: total number of IDE iterations
%
% Written by Heping Song (hepingsong@gmail.com)

path(path, './subfunctions');


%% AÎª¾ØÕó
if ~isa(A, 'function_handle')

Tolerance = 1e-6; % convergence tolerance

%Read the optional inputs
if (rem(length(varargin),2)==1)
    error('Optional inputs must go by pairs!');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case upper('Tolerance')
                Tolerance=varargin{i+1};  
            otherwise
                error(['Unrecognized optional input: ''' varargin{i} '''']);
        end
    end
end


[m, n] = size(A);

if ~exist('opts','var'); opts = []; end

%  options
if isfield(opts,'maxit'); 
    maxit = opts.maxit; 
else
    maxit = 100;
end

if isfield(opts,'x0');   
    x0 = opts.x0;
else
    x0 = zeros(n, 1);
% 	x0 = pinv(A)*y;
end

%%
itr = 1;
x= x0;
r = y;
stop = 0;
AA = pinv(A);


% if isfield(opts,'Th'); 
%     Th = opts.Th; 
% else 
% %     Th = max(abs(x0))/2; 
%     Th = max(abs(AA*y))* beta; 
% end


while ~stop

    % support detection
    actfun = AA * r + x;
    Th = max(abs(actfun))* (beta^itr);
    [support, val] = find(abs(actfun) > Th);
    
    % signal estimation
    % IDE-x approach
    x = zeros(n,1);
    x(support) = pinv(A(:, support)) * y;
    r = y - A*x;    
%     [x r] =MySubsetCG(y, x, A, [], act_set, 1e-16, 0, 100); 
    
    
    
% 	% IDE-s approach
%     x = zeros(n, 1);
%     Aa = A(:, act_set);
%     inact_set = setdiff([1:n], act_set);
%     Ai = A(:, inact_set);
%     P  = pinv(Ai*Ai');
%     xa = pinv(Aa' * P * Aa) * Aa' * y;
%     x(act_set) = xa;
%     x(inact_set) = Ai' * P * (y -  Aa * xa);

    
%     if norm(r) < Tolerance || norm(r)/norm(y) < Tolerance  || itr >100  
	if  norm(r)/norm(y) < Tolerance  || itr >100  
        stop = 1; % convergence tolerance
	else
        itr = itr+1;
%         Th = Th * beta;
	end

end
    
Out.iter=itr;


%%  AÎª¾ä±ú
else   

Tolerance = 1e-10; % convergence tolerance

%Read the optional inputs
if (rem(length(varargin),2)==1)
    error('Optional inputs must go by pairs!');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case upper('At')
                At = varargin{i+1};                
            case upper('MaxIt')
                MaxIt = varargin{i+1};
            case upper('Tolerance')
                Tolerance=varargin{i+1};
%             case upper('Inverse')
%                 AAt_inv = varargin{i+1};
%                 if isa(AAt_inv,'float')
%                    AAt_inv =@(z) AAt_inv*z; 
%                 elseif      ~isa(AAt_inv,'function_handle') 
%                     error('Inverse needs to be a functrion handle or a matrix.')
%                 end
%             case upper('InitialSig')
%                 z_init=varargin{i+1};
%             case upper('Thresh')
%                 thresh=varargin{i+1};
%             case upper('Visibility')
%                 Visibility=varargin{i+1};     
            otherwise
                error(['Unrecognized optional input: ''' varargin{i} '''']);
        end
    end
end


m = length(y);     % measurements length
n = length(At(y)); % signal length


% initialization
t = 1; % iteration number
xt = zeros(n,1);  % initial x
% xt = pinv(A)*y; % initial x
stop = 0; % not convergent
r = y;

AAt_inv = zeros(m,1);
for i=1:m
ee=zeros(m,1);
% ind=randperm(m);
ee(i)=1;
A_row=At(ee);
AAt_inv(i)=1/(norm(A_row)^2);
end

AAt_inv_r = AAt_inv.*r;
zz = xt + At( AAt_inv_r );
z = xt + At( r );
Th = max(abs(zz)) * beta; 

         
while ~stop

%     ee=zeros(m,1);
%     ind=randperm(m);
%     ee(ind(1))=1;
%     A_row=At(ee);
%     AAt_inv=ones(m,1)/(norm(A_row)^2);
%     AAt_inv_y=AAt_inv.*y;

    AAt_inv_r = AAt_inv.*r;
    z = xt + At( AAt_inv_r );
    z = xt + At( r );
    
    % support detection
    [support, val] = find(abs(z) > Th);
    
%     zz = zeros(n,1);
%     zz(act_set) = 1;
    
    % signal estimation
    xt = zeros(n,1);
    [xt r] =MySubsetCG(y, xt, A, At, support, 1e-16, 0, MaxIt); 
    
%     if norm(r) < Tolerance || norm(r)/norm(y) < Tolerance  || t >100  
	if  norm(r)/norm(y) < Tolerance  || t >100  
        stop = 1; % convergence tolerance
    else
        t = t+1;
        Th = Th * beta;
    end
    

end
% fprintf('\n'); 


% output
x=xt;
Out.iter = t;   

  
    
    
    
end
