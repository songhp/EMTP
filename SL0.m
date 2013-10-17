% function s=SL0(A, x, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv, true_s)
function [s  T] =SL0(A, x, varargin)
% SL0(A, x, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv, true_s)
%
%   Returns the sparsest vector s which satisfies underdetermined system of
%   linear equations  A*s=x, using  Smoothed L0  (SL0) algorithm. Note that 
%   the matrix  A should  be a 'wide' matrix  (more columns than rows). The 
%   number of the rows of  matrix A  should  be  equal to the length of the 
%   column vector x.
%
%     The first 3 arguments should necessarily be provided by the user. The 
%   other parameters have defult values  calculated within the function, or
%   may be provided by the user.
%
%   Sequence of Sigma (sigma_min and sigma_decrease_factor):
%     This is a decreasing geometric sequence of positive numbers:
%       - The  first  element   of  the  sequence  of  sigma is  calculated 
%     automatically. The last  element  is  given  by  'sigma_min', and the 
%     change factor for decreasing sigma is given by 'sigma_decrease_factor'. 
%       - The default value of 'sigma_decrease_factor' is 0.5. Larger value 
%     gives better results  for less sparse sources, but it uses more steps 
%     on   sigma   to  reach  sigma_min,  and  hence  it  requires   higher 
%     computational cost.
%       - There is no default  value for  'sigma_min',  and  it  should  be 
%     provided  by  the  user (depending  on his/her estimated source noise 
%     level,  or  his/her  desired  accuracy).  By `noise' we mean here the
%     noise in the sources, that is, the energy of the inactive elements of
%     's'.   For example,  by  the  noiseless  case,  we  mean the inactive
%     elements of 's' are exactly equal to zero. As a rule of tumb, for the
%     noisy case,  sigma_min should be about 2 to 4  times  of the standard
%     deviation of this noise.  For the noiseless case, smaller 'sigma_min'
%     results in  better estimation of the sparsest solution, and hence its
%     value is determined by the desired accuracy.
% 
%   mu_0: 
%        The  value  of  mu_0  scales  the sequence of mu. For each vlue of 
%     sigma, the value of  mu is chosen via mu=mu_0*sigma^2. Note that this 
%     value effects Convergence.
%        The default value is mu_0=2 (see the paper).
%
%   L: 
%        number  of  iterations of the internal (steepest ascent) loop. The
%     default value is L=3.
%
%   A_pinv: 
%        is the  pseudo-inverse of matrix A defined by A_pinv=A'*inv(A*A'). 
%     If it is not provided, it will be calculated within the function.  If
%     you use this function for solving x(t)=A s(t) for different values of
%     't', it would be a good idea to calculate A_pinv outside the function
%     to prevent its re-calculation for each 't'.
%
%   true_s: 
%        is the  true value of the  sparse  solution.  This argument is for
%     simulation purposes. If it is provided by the user, then the function
%     will  calculate the SNR of the estimation for each value of sigma and
%     it provides a progress report.
%
% Authors: Massoud Babaie-Zadeh and Hossein Mohimani
% Version: 1.4
% Last modified: 4 April 2010.
%
%
% Web-page:
% ------------------
%    http://ee.sharif.ir/~SLzero
%
% Code History:
%--------------
% Version 2.0: 4 April 2010
%        Doing a few small modifications that enable the code to work also
%        for complex numbers (not only for real numbers).
%
% Version 1.3: 13 Sep 2008
%        Just a few modifications in the comments
%
% Version 1.2: Adding some more comments in the help section
%
% Version 1.1: 4 August 2008
%    - Using MATLAB's pseudo inverse function to generalize for the case
%      the matrix A is not full-rank.
%
% Version 1.0 (first official version): 4 July 2008.
%
% First non-official version and algorithm development: Summer 2006

% Initialize parameters
sigma_min = 1e-6;
sigma_decrease_factor = 0.5;
mu_0 = 2;
L = 6;

%Read the optional inputs
if (rem(length(varargin),2)==1)
    error('Optional inputs must go by pairs!');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case upper('At')
                At = varargin{i+1};                
            case upper('sigma_min')
                sigma_min = varargin{i+1};
            case upper('sigma_decrease_factor')
                sigma_decrease_factor =varargin{i+1};       
            case upper('mu_0')
                mu_0 = varargin{i+1};                          
            case upper('L')
                L = varargin{i+1};   
            otherwise
                error(['Unrecognized optional input: ''' varargin{i} '''']);
        end
    end
end


% A is a matrix or function handle 
explicitA = ~(ischar(A) || isa(A, 'function_handle'));

if (explicitA)
	A_pinv = pinv(A);
end

%% A is a matrix
if (explicitA)
% Initialization
%s = A\x;
s = A_pinv*x;
sigma = 2*max(abs(s));

T = 0;
% Main Loop
while sigma>sigma_min
    for i=1:L
        delta = OurDelta(s,sigma);
        s = s - mu_0*delta;
        s = s - A_pinv*(A*s-x);   % Projection

    end   
    sigma = sigma * sigma_decrease_factor;
     T = T+1;
end

%% A is a function handle
else
% Initialization
%s = A\x;
s = At(x);
sigma = 2*max(abs(s));

T = 0;
% Main Loop
while sigma>sigma_min
    for i=1:L
        delta = OurDelta(s,sigma);
        s = s - mu_0*delta;
        s = s - At((A(s)-x));   % Projection
    end    
    sigma = sigma * sigma_decrease_factor;
    T = T+1;
end

end %for 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta=OurDelta(s,sigma)

delta = s.*exp(-abs(s).^2/sigma^2);
