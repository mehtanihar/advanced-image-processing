function  res = radon_operator3(m,n,angles1,angles2,angles3)

res.adjoint = 0;
res.row = m;
res.col = n;
res.angles1 = angles1;
res.angles2=angles2;
res.angles3=angles3;

% Register this variable as a partialDCT class
res = class(res,'radon_operator3');
