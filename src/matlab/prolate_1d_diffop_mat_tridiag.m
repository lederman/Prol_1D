function [vdiag, voffdiag] = prolate_1d_diffop_mat_tridiag(c,maxk,isodd)  
%
% Computes the matrix representation of the differential operator,
% in the basis of normalized Legendre polynomials. 
%
% Input:
%   * c : prolate parameters.
%   * maxk  : matrix truncations: the dimensionality of the matrix is k+1
%   * isodd : even or odd submartix
% Output:
%   * vdiag : vector of the elements of the matrix diagonal
%   * voffdiad : vector of coefficients of the matrix off diagonal
%
% Note that the matrix is symmetric tridiagonal. All other elements are
% zeros.
%

%
% Osipov, Andrei, and Vladimir Rokhlin. "On the evaluation of prolate spheroidal wave functions and associated quadrature rules." Applied and Computational Harmonic Analysis 36.1 (2014): 108-142.
%

    numelements = ceil((maxk-isodd+1)/2);
    assert(numelements>1)

    vdiag = zeros(numelements,1);
    voffdiag = zeros(numelements-1,1);
    
    for n=1:numelements
        k = isodd + (n-1)*2;
        vdiag(n,1) = prolate_1d_diffop_mat_diag_element(c,k);
    end
    for n=1:numelements-1
        k = isodd + (n-1)*2;
        voffdiag(n,1) = prolate_1d_diffop_mat_offdiag_element(c,k);
    end
    
end


function v = prolate_1d_diffop_mat_diag_element(c,n)
%
% Osipov, Andrei, and Vladimir Rokhlin. "On the evaluation of prolate spheroidal wave functions and associated quadrature rules." Applied and Computational Harmonic Analysis 36.1 (2014): 108-142.
%

    v = n*(n+1) + (2*n*(n+1)-1)/((2*n+3)*(2*n-1)) *c^2;

end


function v = prolate_1d_diffop_mat_offdiag_element(c,n)
%
% Osipov, Andrei, and Vladimir Rokhlin. "On the evaluation of prolate spheroidal wave functions and associated quadrature rules." Applied and Computational Harmonic Analysis 36.1 (2014): 108-142.
%

    v = ((n+2)*(n+1))/((2*n+3)*sqrt((2*n+1)*(2*n+5))) *c^2;

end

