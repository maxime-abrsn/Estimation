function out = covcorr(data)
% out = covcorr(data)
% data is a matrix of data to calculate the cov/corr matrix
% out is a matrix with covariances below the diag and corr above

% WRITTEN BY: Mike Cliff,  UNC Finance,  mcliff@unc.edu


k = cols(data);
covmat = cov(data);
corrmat = corrcoef(data);

for r = 1:k
  for c = 1:k
    if c > r
      out(r,c) = corrmat(r,c);
    else
      out(r,c) = covmat(r,c);
    end
  end
end

out;
