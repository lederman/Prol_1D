
function [prolate_dat, iserr] = prolate_1d_crea(c,matdim, minEigenvalRatio)

%
% TODO: update fix of eigenvectors, remove dependency on Matlab eig. 
%

    assert(matdim > 4);
    assert(minEigenvalRatio<1)
    assert(minEigenvalRatio>0)

    %
    % Parameters
    %
    iserr = 0;
    prolate_dat.type = 1; 
    prolate_dat.c = c;
    prolate_dat.creaparam.minEigenvalRatio = minEigenvalRatio;
    prolate_dat.creaparam.matdim = matdim;
        
    prolate_dat.evparam.cfs_eps = eps(1.0)/100;

    
    %
    % Evan prolates
    %
    [even_cfs,even_cfs_full,even_lam,even_chi,even_iserr] = prolate_1d_even_crea(c,matdim,minEigenvalRatio);    
    if even_iserr > 0
        iserr = iserr + even_iserr;
    end
    
    %
    % Odd prolates
    %    
    [odd_cfs,odd_cfs_full,odd_lam,odd_chi,odd_iserr] = prolate_1d_odd_crea(c,matdim,minEigenvalRatio);
    if odd_iserr > 0
        iserr = iserr + odd_iserr;
    end
    
    %
    % Combine
    % TODO: add error checking
    %
    raw_chi=zeros(1,max(2*length(even_chi)-1,2*length(even_chi)) );
    raw_chi(1:2:end) = even_chi;
    raw_chi(2:2:end) = odd_chi;
    
    raw_lam=zeros(1,max(2*length(even_lam)-1,2*length(even_lam)) );
    raw_lam(1:2:end) = even_lam;
    raw_lam(2:2:end) = odd_lam;    
    
    raw_cfs_full = zeros( 2*matdim, max(2*length(even_lam)-1,2*length(even_lam)) );
    raw_cfs_full(1:2:end , 1:2:end) = even_cfs;
    raw_cfs_full(2:2:end , 2:2:end) = odd_cfs;
    
    %
    %
    %
    
    ids_to_discard = find( abs(raw_lam)*sqrt(c/(2*pi)) < minEigenvalRatio );
    ids_to_keep = [1:min(ids_to_discard)];
    
    cfs_to_keep = find( max( abs(raw_cfs_full) , [], 2) > prolate_dat.evparam.cfs_eps );
    cfs_to_keep = [1:max(cfs_to_keep)];
        
    %
    % save
    %
    prolate_dat.chi = raw_chi(ids_to_keep);
    prolate_dat.lam = raw_lam(ids_to_keep);
    prolate_dat.nu  = prolate_dat.lam *sqrt(c/(2*pi));
    
    prolate_dat.cfs = raw_cfs_full( cfs_to_keep, ids_to_keep );
    
    
    prolate_dat.even.cfs = prolate_dat.cfs(1:2:end , 1:2:end);
    prolate_dat.odd.cfs = prolate_dat.cfs(2:2:end , 2:2:end);
    
    prolate_dat.num_prols = ids_to_keep(end);
    
    %
    % Warnings
    % Take plenty of margin for the truncation.
    %
    if (ids_to_keep(end)+20 >= 2*matdim)
        warning('prolate_1d_crea: insufficient margin in matrix size (number of prolates)')
        iserr = iserr+10;
    end
    if (cfs_to_keep(end)+40 >= 2*matdim)
        warning('prolate_1d_crea: insufficient margin in matrix size (number of coefficients)')
        iserr = iserr+100;        
    end

    
end



function [cfs,cfs_full,lam,chi,iserr] = prolate_1d_odd_crea(c,matdim,minEigenvalRatio)
    
    assert(matdim>1)
    
    iserr = 0;
    isodd = 1;
    maxk = 2*(matdim-1)+1;
    [mat, vdiag, voffdiag ] = prolate_diffop_mat_full(c,maxk,isodd);

    % TODO: remove dependency on matlab eig!
    [u,d] = eig(mat);
    % Note that matlab keeps the first small elements in the eigenvector,
    % but this is an undocumented feature.
    % Note that Matlab eig truncates some small coefficients by setting them to 0. 
    [chi,eigvals_order] = sort(diag(d),'ascend');
    chi=chi.';
    cfs = u(:,eigvals_order);
    cfs_full = zeros( 2*size(cfs,1) , size(cfs,2) );
    cfs_full(2:2:end,:) = cfs;
    
    %This is the code if there were no truncation in matlab eig:    
    %[prol_at_0, dprol_at_0] = prolate_LegeNorm_ex(cfs_full, 0.0 ) ;    
    %lam = 1i * sqrt(2/3) * c * cfs(1,:) ./ dprol_at_0;
    % % standard sign:
    %cfs = bsxfun(@times, cfs, sign(dprol_at_0));
    %cfs_full = bsxfun(@times, cfs_full, sign(dprol_at_0));
    
    
    %
    % fix each vector, due to matlab truncation
    % TODO: replace with independent implementation
    %
    
    %raw_cfs = cfs; % for debugging
    lam = zeros(1,matdim);
    fix_niter = 5; 
    for jn=1:matdim-1
        cfs(:,jn) = prolate_1d_crea_fix_eigenvec(mat, cfs(:,jn) , chi(jn), chi(jn+1) , fix_niter);
        cfs_full(1:2:end,jn) = cfs(:,jn);
        [prol_at_0,dprol_at_0] = prolate_LegeNorm_ex(cfs_full(:,jn), 0.0 ) ;   
        lam(jn) =  1i * sqrt(2/3) * c * cfs(1,jn) ./ dprol_at_0;
        cfs(:,jn) = cfs(:,jn)* sign(dprol_at_0);
        cfs_full(:,jn) = cfs_full(:,jn) * sign(dprol_at_0);
        if ( abs(lam(jn))*sqrt(c/(2*pi)) < minEigenvalRatio )
            break;
        end
    end
    if (jn+20 > matdim)
        iserr = 10;
        warning('prolate_1d_odd_crea : matdim probably too small');
    end
    
    cfs = cfs(:,1:jn);
    cfs_full = cfs_full(:,1:jn);
    lam = lam(:,1:jn);
    chi = chi(:,1:jn);
    
    
end



function [cfs,cfs_full,lam,chi,iserr] = prolate_1d_even_crea(c,matdim,minEigenvalRatio)
    
    assert(matdim>1)
    
    iserr = 0;
    isodd = 0;
    maxk = 2*(matdim-1);
    [mat, ~, ~ ] = prolate_diffop_mat_full(c,maxk,isodd);

    % TODO: remove dependency on matlab eig!
    [u,d] = eig(mat);
    % Note that matlab keeps the first small elements in the eigenvector,
    % but this is an undocumented feature.
    % Note that Matlab eig truncates some small coefficients by setting them to 0. 
    [chi,eigvals_order] = sort(diag(d),'ascend');
    chi=chi.';
    cfs = u(:,eigvals_order);    
    cfs_full = zeros( 2*size(cfs,1)-1 , size(cfs,2) );
    cfs_full(1:2:end,:) = cfs;
    
    %This is the code if there were no truncation in matlab eig:
    %prol_at_0 = prolate_LegeNorm_ex(cfs_full, 0.0 ) ;    
    %lam = sqrt(2) * cfs(1,:) ./ prol_at_0;
    %cfs = bsxfun(@times, cfs, sign(prol_at_0));
    %cfs_full = bsxfun(@times, cfs_full, sign(prol_at_0));
    
    %
    % fix each vector, due to matlab truncation
    % TODO: replace with independent implementation
    %
    
    %raw_cfs = cfs; % for debugging
    lam = zeros(1,matdim);
    fix_niter = 5;
    for jn=1:matdim-1
        cfs(:,jn) = prolate_1d_crea_fix_eigenvec(mat, cfs(:,jn) , chi(jn), chi(jn+1),fix_niter);
        cfs_full(1:2:end,jn) = cfs(:,jn);
        prol_at_0 = prolate_LegeNorm_ex(cfs_full(:,jn), 0.0 ) ;   
        lam(jn) = sqrt(2) * cfs(1,jn) ./ prol_at_0;
        cfs(:,jn) = cfs(:,jn)* sign(prol_at_0);
        cfs_full(:,jn) = cfs_full(:,jn) * sign(prol_at_0);
        if ( abs(lam(jn))*sqrt(c/(2*pi)) < minEigenvalRatio )
            break;
        end
    end
    if (jn+20 > matdim)
        iserr = 10;
        warning('prolate_1d_even_crea : matdim probably too small');
    end
    
    cfs = cfs(:,1:jn);
    cfs_full = cfs_full(:,1:jn);
    lam = lam(:,1:jn);
    chi = chi(:,1:jn);
    
    
end



function v = prolate_1d_crea_fix_eigenvec(mat, vec0 , chi, chinext, fix_niter)
    
    tmpmat = mat-eye(size(mat,1))*...
        (chi+(chinext-chi)/10^5 );
    v=vec0;
    for jn=1:fix_niter
        v=tmpmat\v;
        v = v/norm(v);
    end
    
end

