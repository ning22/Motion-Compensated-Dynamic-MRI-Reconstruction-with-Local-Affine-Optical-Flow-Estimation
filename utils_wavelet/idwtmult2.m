% idwtmult2.m
%
% Inverts dwtmult2.
% Usage : x = idwtmult2(w, g0, g1, L, sym)
%
% Written by : Justin Romberg
% Created : 6/26/2001

function x = idwtmult2(w, g0, g1, L, sym)

if (nargin == 4), sym = 0; end

[Nr,Nc] = size(w);

for ll = L:-1:1
  for jj = 1:(Nc*2^(-ll+1))
    w(1:Nr*2^(-ll+1),jj) = idwtlevel1(w(1:Nr*2^(-ll+1),jj)', g0, g1, sym)';
  end
  for ii = 1:(Nr*2^(-ll+1))
    w(ii,1:Nc*2^(-ll+1)) = idwtlevel1(w(ii,1:Nc*2^(-ll+1)), g0, g1, sym);
  end
end
x = w;

