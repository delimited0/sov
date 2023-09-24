n = 4;
a = -2*ones(n, 1);
b = ones(n, 1);
cn = eye(n);

Sigma = .5 * eye(n) + .5 * ones(n, 1) * ones(1, n);

[V, D] = eig(Sigma);
[d,pm] = sort(diag(D)','descend');
d = sqrt(max(d,0));  
np = sum(d>0); 
c = V(:,pm(1:np)).*( ones(n,1)*d(1:np) ); 
ch = ch*c; 



r = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5];
a = -inf(4,1); b = [ 1 2 3 4 ]';
p = mvnxpb( r, a, b ); disp(p)

p = bvnu( -3, -1, .35 )

p = bvnu(


%% lattice rule comparison
[lat, n] = fstrnk(99991, 100);