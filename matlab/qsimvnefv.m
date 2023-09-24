function [ p, e, ef, efe ] = qsimvnefv( m, r, a, b, f )
% 
%    Qsimvnefv is a vectorized version of qsimvnef which uses a
%    randomized quasi-random rule with m points to estimate an MVN 
%    expectation for a positive definite covariance matrix r, with lower
%    integration limits a and upper integration limits b, and expectation 
%    function f. F, which may be a (column) vector function, must be 
%    vectorized, with a matrix input having f input components as columns.  
%   Probability p is output with error estimate e, along with
%      expected value(s) ef and error estimate(s) efe.
%    Example use:
%     r = [ 4 3 2 1; 3 5 -1 1; 2 -1 4 2; 1 1 2 5 ];
%     a = -inf*[ 1 1 1 1 ]'; b = [ 1 2 3 4 ]';
%%   To compute expected values for all four variables and their squares
%     f = @(x)[ x(1:4,:); x(1:4,:).^2 ]; 
%        note required vector structure for x and f components.
%     [ p e ef efe ] = qsimvnefv( 50000, r, a, b, f ); disp([[p e];[ef efe]])
%   Note: if f is defined in an m-file f.m, then use
%     [p e ef efe] = qsimvnefv( 50000, r, a, b, 'f' ); disp([[p e];[ef efe]])
%
%   This function uses an algorithm given in the paper
%      "Numerical Computation of Multivariate Normal Probabilities", in
%      J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
%          Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
%          Email : AlanGenz@wsu.edu
%  The primary references for the numerical integration are 
%   "On a Number-Theoretical Integration Method"
%   H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11, and
%   "Randomization of Number Theoretic Methods for Multiple Integration"
%    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13(1976), pp. 904-14.
%
%   Alan Genz is the author of this function and following Matlab functions.
%
%
% Copyright (C) 2013, Alan Genz,  All rights reserved.               
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided the following conditions are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the 
%      distribution.
%   3. The contributor name(s) may not be used to endorse or promote 
%      products derived from this software without specific prior 
%      written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
% OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% Initialization
%
[n, n] = size(r); ch = chol(r)'; ct = ch(1,1); ai = a(1); bi = b(1); 
ci = ( 1 + sign(ai) )/2; if abs(ai) < 9*ct, ci = Phi(ai/ct); end,  
di = ( 1 + sign(bi) )/2; if abs(bi) < 9*ct, di = Phi(bi/ct); end, di = di - ci; 
ps = sqrt( primes( 5*(n+1)*log(n+1)/4 ) ); q = ps(1:n)'; % Richtmyer generators
ns = 10; nv = max( m/ns, 1 ); p = 0; e = 0; ef = 0; efe = 0; on = ones(1,nv);
lf = length( f((max(a,-9)+min(b,9))/2) ); y = zeros(n,nv); 
%
% Randomization loop for ns samples
%
for j = 1 : ns, c = ci*on; dc = di*on; pv = dc; 
  for i = 2 : n, x = abs( 2*mod( q(i-1)*[1:nv] + rand, 1 ) - 1 ); 
    y(i-1,:) = Phinv( c + x.*dc ); s = ch(i,1:i-1)*y(1:i-1,:); 
    ct = ch(i,i); ai = a(i) - s; bi = b(i) - s; c = on; d = c;
    c(find( ai < -9*ct )) = 0; d(find( bi < -9*ct )) = 0; 
    tstl = find( abs(ai) < 9*ct ); c(tstl) = Phi(ai(tstl)/ct);
    tstl = find( abs(bi) < 9*ct ); d(tstl) = Phi(bi(tstl)/ct);
    dc = d - c; pv = pv.*dc; 
  end, pj = mean(pv); y(n,:) = Phinv(c+abs(2*mod(q(n)*[1:nv]+rand,1)-1).*dc);  
  y = min(max(-9,y),9); efj = mean( repmat(pv,lf,1).*f(ch*y), 2 ); 
  dp = (  pj - p  )/j;  p =  p + dp;   e =   e*( j - 2 )/j +  dp^2; 
  df = ( efj - ef )/j; ef = ef + df; efe = efe*( j - 2 )/j + df.^2;
end, e = 3*sqrt(e); efe = 3*sqrt(efe); % error estimates = 3 x standard errors
if p > 0, ef = ef/p; efe = efe/p; end
%
% end qsimvnefv
%
%
%  Standard statistical normal distribution functions
%
function p =   Phi(z), p =  erfc( -z/sqrt(2) )/2;
% function z = Phinv(p), z = -sqrt(2)*erfcinv( 2*p );
function z = Phinv(p), z = norminv(p); % for newer Matlabs, or Octave
%


