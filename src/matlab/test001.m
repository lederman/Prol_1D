function test001()


    
    c=50;
    matdim=200;
    minEigenvalRatio = 10^-40;
    
    [prolate_dat, iserr] = prolate_1d_crea(c,matdim, minEigenvalRatio);

    figure; 
    semilogy((abs(prolate_dat.nu)))
    ylim([10^-30,3])
    title('magnitude of eigenvalues (not scaled)')
    
    xx=linspace(-1,1,500);
    prolate_ids = [1:10];
    [v] = prolate_1d_ev(prolate_dat, prolate_ids, xx);
    
    figure; 
    plot(v(:,1:2:end));
    title('even prolates')
    figure; 
    plot(v(:,2:2:end));
    title('odd prolates')
    
    
   
  
end







