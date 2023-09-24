function [ z, n ] = fstrnk( ni, sm, om, gm, bt )
% 
% Reference: 
%   "Fast Component-by-Component Construction, a Reprise for Different 
%     Kernels", D. Nuyens and R. Cools. In H. Niederreiter and D. Talay,
%     editors, Monte-Carlo and Quasi-Monte Carlo Methods 2004, 
%     Springer-Verlag, 2006, 371-385.
% Modifications to original by A. Genz, 05/07
% Typical Use:  
%  om = inline('x.^2-x+1/6'); n = 99991; s = 100; gam = 0.9.^[1:s];
%  z = fastrank( n, s, om, gam, 1 + gam/3 ); disp([z'; e])
%
n = fix(ni); if ~isprime(n), pt = primes(n); n = pt(length(pt)); end
if nargin < 3, om = inline('x.^2-x+1/6'); 
  bt = ones(1,sm); gm = [ 1 (4/5).^[0:sm-2] ]; 
end, q = 1; w = 1; z = [1:sm]'; m = ( n - 1 )/2; g = prmrot(n); 
perm = [1:m]; for j = 1 : m-1, perm(j+1) = mod( g*perm(j), n ); end
perm = min( n - perm, perm ); c = om(perm/n); fc = fft(c); 
for s = 2 : sm, q = q.*( bt(s-1) + gm(s-1)*c([w:-1:1 m:-1:w+1]) );
  [ es w ] = min( real( ifft( fc.*fft(q) ) ) ); z(s) = perm(w); 
end