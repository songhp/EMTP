function [xt out] = EMTP(A, y, s, mu, varargin)
% EMTP: Expectation conditinal Maximazation either Thresholding Pursuit
% Coded by Heping Song (hepingsong@gmail.com)
% 
%%%%%%%%%%%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT(compulsory):
%   A: m*n measurement matrix
%   y: m*1 measurements
%   s: step size for estimeate the sparsity, recommend m/(4*log(n))
%   
% INPUT(optional):
% 'mu': normalized term for x + mu * AA * r
% 'At': the adjoint function handle that computes A^T*y for function handle

%%%%%%%%%%%%%%%%%% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   xt: reconstructed sparse signal
%   out.iter: number of iterations

path(path, './subfunctions');


%% AÎª¾ØÕó
if ~isa(A, 'function_handle')

Tolerance = 1e-10; % convergence tolerance

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


n = size(A,2); % signal length

% initialization
t = 1; % iteration number
xt = zeros(n,1);  % initial x
% xt = pinv(A)*y; % initial x
r = y;  % initial residue
support_size = s;  % initial sparsity level
support = [];
stop = 0; % not convergent

AA = pinv(A);

while ~stop

    % support detection
    [val idx] = sort(abs(xt + mu * AA *r), 'descend');
	support = sort(idx(1:support_size));
    
    % signal estimation
    xt = zeros(n,1);
    xt(support) = pinv(A(:, support)) * y;
    r = y - A*xt;
    
%     if norm(r) < Tolerance || norm(r)/norm(y) < Tolerance  || t >100  
	if  norm(r)/norm(y) < Tolerance  || t >100 
        stop = 1; % convergent
    else
        t = t+1;
    end
    

end
% fprintf('\n'); 
% % reconstruction
% xt = zeros(n,1);
% xt_support = pinv(A(:,support))*y;
% xt(support) = xt_support;


% output
out.iter = t;

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
r = y;  % initial residue
support_size = s;  % initial sparsity level
support = [];
stop = 0; % not convergent


ee=zeros(m,1);
ind=randperm(m);
ee(ind(1))=1;
A_row=At(ee);
AAt_inv=ones(m,1)/(norm(A_row)^2);
AAt_inv_y=AAt_inv.*y;


while ~stop

    Ax_hat =A(xt);
    AAt_invAx_hat =AAt_inv.*Ax_hat;

    z = xt + At( AAt_inv_y - mu* AAt_invAx_hat );
    % support detection
    [val idx] = sort(abs(z), 'descend');
	support = sort(idx(1:support_size));
    zz = zeros(n,1);
    zz(support) = 1;
    
    % signal estimation
    xt = zeros(n,1);
    [xt r] =MySubsetCG(y, xt, A, At, find(zz~=0), 1e-9, 0, MaxIt); 
    
%     if norm(r) < Tolerance || norm(r)/norm(y) < Tolerance  || t >100
    if  norm(r)/norm(y) < Tolerance  || t >100
        stop = 1; % convergence tolerance
    else
        t = t+1;
    end
    

end
% fprintf('\n'); 
% % reconstruction
% xt = zeros(n,1);
% xt_support = pinv(A(:,support))*y;
% xt(support) = xt_support;


% output
out.iter = t;   
    
end    

