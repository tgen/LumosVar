function y = bbinopdf_ln(x,n,a,b)
%BBINOPDF Beta-binomial probability distribution function.
% Y = BBINOPDF_LN(X,N,A,B) returns the beta-binomial probability
% density function with parameters N, A and B at the values in X. 
% Note: The density function is zero unless N, A and B are integers.
%
% The Beta-binomial distribution is used to model the number of successes
% in n binomial trials when the probability of success p is a Beta(a,b)
% random variable. The extreme flexibility of the shape of the Beta 
% distribution means that it is often a very fair representation of the 
% randomness of p.
%
% A variable with a Beta-binomial distribution is distributed as a binomial
% distribution with parameter p, where p is distribution with a beta 
% distribution with parameters a (alpha) and b (beta). For n trials, it has
% probability density function:
%
%                                 B(x+a,n-x+b)
%                   p(x) = n_C_x ---------------
%                                     B(a,b)
%
% where B(a,b) is a beta function and n_C_x is a binomial coefficient.
% (http://mathworld.wolfram.com/BetaBinomialDistribution.html)
% (http://en.wikipedia.org/wiki/Beta-binomial_model)
%
% The probability of success varies randomly, but in any one scenario that
% probability applies to all trials. For example, you might consider using
% the Beta-binomial distribution to model:
%
% --The number of life insurance policy holders who will die in any one 
% year, where some external variable (e.g. highly contagious disease, 
% extreme weather) moderates the probability of death of all individual to
% some degree.
% --The number of cars that crash in a race of n cars, where the predominant
% factor is not the skill of the individual driver, but the weather on the
% day.
% --The number of bottles of wine from a producer that are bad where the 
% predominant factor is not how each bottle is treated, but something to do
% with the batch as a whole.
%
% The Beta-binomial is a two-dimensional multivariate Polya distribution, 
% as the binomial and beta distributions are special cases of the 
% multinomial and Dirichlet distributions, respectively.
%
% Syntax: function y = bbinopdf_ln(x,n,a,b) 
%      
% Inputs:
% x - number of success  
% n - number of trials  
% a - alpha parameter of Beta
% b - beta parameter of Beta
%
% Output:
% y - beta-binomial probability value
%
% Example. Suppose we make 5 trials on sample with a Beta-binomial 
% distribution with parameters alpha = 3, and beta = 4, and we are 
% interested to get the probability of get exactly 3 successes.
%
% Calling on Matlab the function: 
%             y = bbinopdf_ln(3,5,3,4)
%
% Answer is: (in format long)
%
% y =
%
%   0.21645021645022
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, F.A. Trujillo-Perez
%            and N. Castro-Castro
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
% Copyright. September 28, 2009.
%
% ---Con cariño para Ney.----
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, F.A. Trujillo-Perez and 
%   N. Castro-Castro (2009). bbinopdf:Beta-binomial probability 
%   densiy function. A MATLAB file. [WWW document]. URL 
%   http://www.mathworks.com/matlabcentral/fileexchange/25454
%
% Modified by Rebecca Halperin to use log transformation to avoid numeric overflow
% Last Modified June 8th 2016

if nargin < 4, 
    error('bbinopdf:TooFewInputs','Requires four input arguments.'); 
end
% 
% [errorcode x n a b] = distchck(4,x,n,a,b);
% 
% if errorcode > 0
%     error('bbinopdf:InputSizeMismatch',...
%         'Requires non-scalar arguments to match in size.');
% end
% 
% if (length(x)~=1) || (fix(x) ~= x) || (x < 0),
%    error('bbinopdf:InvalidData',...
%        'BBINOPDF requires that X must be a non-negative and integer.')
% end
% 
% if (length(n)~=1) || (fix(n) ~= n) || (n < 0),
%    error('bbinopdf:InvalidData',...
%        'BBINOPDF requires that N must be a non-negative and integer.')
% end
% 
% if (length(a)~=1) || (fix(a) ~= a) || (a < 0),
%    error('bbinopdf:InvalidData',...
%        'BBINOPDF requires that A must be a non-negative and integer.')
% end
% 
% if (length(b)~=1) || (fix(b) ~= b) || (b < 0),
%    error('bbinopdf:InvalidData',...
%        'BBINOPDF requires that B must be a non-negative and integer.')
% end

y = exp(gammaln(n + 1)-gammaln(x + 1)-gammaln(n - x + 1)+betaln((a + x),(b + n - x))-betaln(a,b));

%alternatively, using the Gamma function, it can be expressed as:
%
%y = exp(gammaln(n + 1)-gammaln(x + 1)-gammaln(n - x + 1))*(gamma(a + b)*...
%    gamma(a + x)*gamma(b + n - x))/(gamma(a)*gamma(b)*gamma(a + b + x))

return,
