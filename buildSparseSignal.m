function x = buildSparseSignal(dimensions,sparsity,distributiontype,seed)

randn('seed',seed);
rand('seed',seed);
x = zeros(dimensions,1);
activelements = randperm(dimensions);

switch lower(distributiontype)
	case 'normal'
		v = randn(sparsity,1);
	case 'bimodalgaussian'
		v = sign(randn(sparsity,1)).*(3+randn(sparsity,1));
	case 'laplacian'
		v = randlap([sparsity,1],10,seed);
	case 'uniform'
		v = sign(randn(sparsity,1)).*rand(sparsity,1);
	case 'bimodaluniform'
		v = sign(randn(sparsity,1)).*(3+rand(sparsity,1));
	case 'bernoulli'
		v = sign(randn(sparsity,1));
    case 'bimodalrayleigh'
    v = 3*sign(randn(sparsity,1)).*sqrt(-2*log(1-rand(sparsity,1)));
    case 'powerlaw'
        v = sign(randn(sparsity,1)).*((1:sparsity).^(-1/0.5))'; % lambda=0.5
    case 'exponential'
        v = sign(randn(sparsity,1)).*(exp(-0.5*(1:sparsity)))'; % lambda=0.5
end
x(activelements(1:sparsity)) = v;


function x = randlap(siz,lambda,seed)

rand('seed',seed);

% RANDL  random numbers distributed according to the Laplace distribution
%   RANDL(N) will return an NxN matrix containing pseudo-random values
%   drawn from a Laplace distribution with zero mean and standard deviation
%   one. RAND([M,N,...,P]) will return an MxNx...xP matrix.
%   RANDL([M,N,...,P],lambda) will return an MxNx...xP matrix of
%   pseudo-random numbers with parameter lambda. CAUTION: the pdf is
%   assumed as 
%           pdf = lambda/2 * exp(-lambda*abs(x))
%
% The Laplace random numbers are generated using the the RAND function to
% generate uniformly distributed numbers in (0,1) and then the probability
% integral transformation, i.e. using the fact that 
%   if Fx(x) is the cdf of a random variable x, then the RV z=Fx(x) is
%   uniformly distributed in (0,1).
%
% In order to generate random numbers with mean mu and variance v, then
% generate the random numbers X with mean 0 and variance 1 and then
%       X = X*sqrt(v)+mu;
%
% C. Saragiotis, Oct 2008


if nargin==1, 
    lambda = sqrt(2); % this gives a std=var=1.
end

z = rand(siz);
x = zeros(siz);

in = z<=.5;
ip = z> .5;
x(in) =  1/lambda *log(2*z(in));
x(ip) = -1/lambda *log(2*(1-z(ip)));