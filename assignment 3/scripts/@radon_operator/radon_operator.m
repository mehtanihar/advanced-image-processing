function  res = radon_operator(m,n,angles)

res.adjoint = 0;
res.row = m;
res.col = n;
res.angles = angles;

% Register this variable as a partialDCT class
res = class(res,'radon_operator');
