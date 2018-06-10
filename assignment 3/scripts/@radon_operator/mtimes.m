function res = mtimes(A,x)

if A.adjoint == 0 %A*x
    res = reshape((radon(dct2(reshape(x,[A.row,A.col])),A.angles)),[],1);%MN->LP
    
else %At*x
    res = reshape((idct2(iradon(reshape(x,[],size(A.angles,1)),A.angles,'linear','Ram-Lak',1,A.row))),A.row*A.col,1);
end
