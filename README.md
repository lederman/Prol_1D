# Prol 1D
Prolate Spheroidal Wave Functions (PSWF)

In development.

Matlab version. 
This implementation does not scale well at this point due to inefficient use of Matlab's eigendecomposition, which will be replaced in future versions: 
* Matlab's eig takes full matrices (although it exploits the tridiagonal structure internally).
* Matlab truncates certain small elements of the eigenvectors. Matlab is accurate - this truncation is well within the specs of eigendecomposition, in the sense that the vectors are accurate in l^2 norm, but alternative algorithms can retain these elements, which are very useful in prolate numerical algorithms. Currently, we bypass this problem using a slow inverse power method. 
* There are ways to exploit analytical results to make the eigendecomposition faster. 

