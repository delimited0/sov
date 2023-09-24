function r = prmrot(pin)
%
% find primitive root for prime p, might fail for large primes (p > 32e7)
%
p = pin; if ~isprime(p), pt = primes(p); p = pt(length(pt)); end
pm = p - 1; fp = unique(factor(pm)); n = length(fp); r = 2; k = 1;
while k <= n; d = pm/fp(k); rd = r;
  for i = 2 : d, rd = mod( rd*r, p ); end % computes r^d mod p
  if rd == 1, r = r + 1; k = 0; end, k = k + 1;
end 
