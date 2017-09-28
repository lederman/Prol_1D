
function [v ,dv] = prolate_Lege_ex(cfsvec,xx)   
%
% Evaluates functions and their derivatives expanded in Legendre Polynomials  P_n(x) (not normalized)
%
%  v(i,j) = \sum_{q=0}^{k-1} cfsvec(q+1,j) P_q(x_i)
%
% Input:
%   * cfsvec : k x m matrix. 
%       Columns of coefficients, each column has the k coefficients of an expansion 
%       in Legendre polynoimals of order k-1 for one of the m different functions.
%   * xx :  a vector of length l. 
%       each entry is a value of x where each one of the k expansions should be evaluated.
% Output:
%   * v : l x m matrix.
%       The j-th column is the j-th function evaluated at the l points.
%   * dv : derivative
%

    xx = xx(:); 

    jcm1 = 0*xx(:)+1;    
    jc0  = xx(:);
    jcp1 = 0*xx(:);    
    if (size(cfsvec,1)>0)
        v= jcm1 *cfsvec(1,:);         
    else
        v=zeros([length(xx),size(cfsvec,2)],class(xx));
    end
    dv=0*v;
    djcm1  = 0*xx(:);    
    djc0  = 0*xx(:)+1;    
    djcp1  = 0*xx(:); 
    if size(cfsvec,1)>1        
        v=v+jc0 * cfsvec(2,:);
        dv=dv + djc0 * cfsvec(2,:);
        djcp1(:) = djcm1(:) + (2*1+1)*jc0(:);
        djcm1(:) = djc0;
        djc0(:)=djcp1(:);
        for k=1:size(cfsvec,1)-2
            jcp1(:) = (((2*k+1)*xx(:)).*jc0 -(k)*jcm1) / ((k+1));
            jcm1(:) = jc0(:);
            jc0(:)  = jcp1(:);            
            v=v+jc0*cfsvec(k+2,:);
            dv=dv + djc0 * cfsvec(k+2,:);
            djcp1(:) = djcm1(:) + (2*(k+1)+1)*jc0(:);
            djcm1(:) = djc0;
            djc0(:)=djcp1(:);
        end
    end
end
