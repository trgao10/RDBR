% prepare data
Xi = (1:5)';
X = (0:0.1:6)';
sigma = 1;
W = exp(-getEuclideanDistance(X',Xi').^2/sigma^2);
ord = 2;
tol = 1e-8;
% get MLS coefficients
coeff = MLS(Xi,X,W,ord,tol);
% randomly sample points
fi = rand(size(Xi));
% interpolate
f = coeff*fi;
% plot
figure;
plot(Xi,fi,'r*');
hold all
plot(X,f,'b');
