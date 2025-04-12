function J = comp_current(r,mu,Vth,z,n, Bp, Bn)
  % This computes the current xj(x) at the cell centers 
  % hence it is of size lr-1
  % can be used to compute the total current 
  % or be integratted with the constant piecewise method

  nBn = n(2:end).*Bn;
  nBp = n(1:end-1).*Bp;

  J = z*mu*Vth./(log(r(1:end-1) ./ r(2:end))).*(nBn-nBp);
end

    










