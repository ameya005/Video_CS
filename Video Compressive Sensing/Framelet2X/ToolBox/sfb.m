function y = sfb(lo, hi1, hi2, sf)

% Synthesis filter bank
%
% USAGE:
%    y = sfb(lo, hi1, hi2, sf)
% INPUT:
%    lo - lowpass input
%    hi1 - bandpass input
%    hi2 - highpass input
%    sf - synthesis filters
% OUTPUT:
%    y - output signal
% See also afb
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

N = 2*length(lo);

g0 = sf{1};
g1 = sf{2};
g2 = sf{3};

L0 = length(g0);
L1 = length(g1);
L2 = length(g2);

% symmetric extension
A = L0/2;
lo = [lo(A:-1:1) lo lo(N/2:-1:N/2-A)];
% lowpass filter
lo = up(lo, 2);
lo = conv(lo, g0);
% extract valid part
lo = lo(3*L0/2-1+[1:N]);

% normalize
hi1(1) = sqrt(2)*hi1(1);
hi1(N/2+1) = sqrt(2)*hi1(N/2+1);
% symmetric extension
A = L1/2;
hi1 = [hi1(A:-1:2) hi1 hi1(N/2:-1:N/2-A)];
% highpass filter
hi1 = up(hi1, 2);
hi1 = conv(hi1, g1);
% extract valid part
hi1 = hi1(3*L1/2-2+[1:N]);

% symmetric extension
A = L2/2;
hi2 = [0 hi2 0];
hi2 = [-hi2(A:-1:2) hi2 -hi2(N/2:-1:N/2-A)];
% highpass filter
hi2 = up(hi2, 2);
hi2 = conv(hi2, g2);
% extract valid part
hi2 = hi2(3*L1/2-2+[1:N]);

% add signals
y = lo + hi1 + hi2;
