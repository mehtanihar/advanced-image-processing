function  res = radon_operator2(m,n,angles1,angles2)

res.adjoint = 0;
res.row = m;
res.col = n;
res.angles1 = angles1;
res.angles2=angles2;

% Register this variable as a partialDCT class
res = class(res,'radon_operator2');
