function res = mtimes(A,x)

if A.adjoint == 0 %A*x
    A1=radon_operator(A.row,A.col,A.angles1);
    A2=radon_operator(A.row,A.col,A.angles2);
    res=[A1*x(1:A.row*A.col,1);A2*(x(1:A.row*A.col,1)+x(A.row*A.col+1:2*A.row*A.col,1))];%2MN->2LP
    
else %At*x
    A1=radon_operator(A.row,A.col,A.angles1);
    A2=radon_operator(A.row,A.col,A.angles2);
    res=[A1'*x(1:size(x,1)/2,1);A2'*x(size(x,1)/2+1:size(x,1),1)-A1'*(x(1:size(x,1)/2,1))];%2LP->2MN
    
end
