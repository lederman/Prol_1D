
function [mat, vdiag, voffdiag ] = prolate_diffop_mat_full(c,maxk,isodd)
%
%


    [vdiag, voffdiag] = prolate_1d_diffop_mat_tridiag(c,maxk,isodd);
    mat = zeros(length(vdiag));
    
    for k=1:length(vdiag)
        mat(k,k) = vdiag(k);        
    end
    for k=1:length(voffdiag);
        mat(k,k+1) = voffdiag(k);        
        mat(k+1,k) = voffdiag(k);        
    end
    
end

