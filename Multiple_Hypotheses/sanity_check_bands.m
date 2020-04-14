function sanity_check_bands(Pmin,Pmax,dx)

if any(size(Pmin) ~= size(Pmax))
    error('Densities must be of the same size.');
elseif any(sum(Pmin,2)*dx > 1) || any(sum(Pmax,2)*dx < 1) || any(Pmin(:) < 0) || any(Pmin(:) > Pmax(:))
    error('Invalid corridor specifiactions.');
end
