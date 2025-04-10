function J = comp_current(r,mu,v,Vth,z,n)
  % This computes the current xj(x) at the cell centers 
  % hence it is of size lr-1
  % can be used to compute the total current 
  % or be integratted with the constant piecewise method

  DV = z*diff(v) / Vth;
  [Bp, Bn] = bimu_bernoulli (DV);

  nBn = n(2:end).*Bn;
  nBp = n(1:end-1).*Bp;

  % unsure of the z
  J = z*mu*Vth./(log(r(1:end-1) ./ r(2:end))).*(nBn-nBp);
end

    










