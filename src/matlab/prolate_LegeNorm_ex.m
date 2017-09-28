
function [v ,dv] =  prolate_LegeNorm_ex(cfsvec,xx) 
%
% Evaluates functions and their derivative expanded in NORMALIZED Legendre Polynomials  \hat{P}_n(x) 
%
%  v(i,j) = \sum_{q=0}^{k-1} cfsvec(q+1,j) \hat{P}_q(x_i)
%
% Input:
%   * cfsvec : k x m matrix. 
%       Columns of coefficients, each column has the k coefficients of an expansion 
%       in normalized Legendre polynoimals of order k-1 for one of the m different functions.
%   * xx :  a vector of length l. 
%       each entry is a value of x where each one of the k expansions should be evaluated.
% Output:
%   * v : l x m matrix.
%       The j-th column is the j-th function evaluated at the l points.
%   * dv : derivative
%
%    
   
    K=size(cfsvec,1)-1;
    cfsvec_lege = cfsvec;
    cfsvec_lege = bsxfun(@times, sqrt(1/2*(2*[0:K]' + 1)), cfsvec_lege);
    
    [v ,dv] = prolate_Lege_ex(cfsvec_lege,xx) ;
    
 
end


