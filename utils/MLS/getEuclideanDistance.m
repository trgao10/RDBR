function r = getEuclideanDistance(u,v)

r = sqrt(bsxfun(@plus,sum(u.^2,1)',sum(v.^2,1)) - 2*u'*v);
