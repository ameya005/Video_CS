function [lo, hi1, hi2] = afb(x, af)

% Analysis filter bank
%
% USAGE:
%    [lo, hi1, hi2] = afb(x, af)
% INPUT:
%    x - N-point vector, with N even
%    af - analysis filters
%    af{1} - lowpass filter (symmetric even-length)
%    af{2} - bandpass filter (symmetric even-length)
%    af{3} - highpass filter (antisymmetric even-length)
% OUTPUT:
%    lo - lowpass output
%    hi1 - bandpass output
%    hi2 - highpass output
% EXAMPLE:
%    [af, sf] = filters1;
%    x = rand(1,64);
%    [lo, hi1, hi2] = afb(x, af);
%    y = sfb(lo, hi1, hi2, sf);
%    err = x - y; 
%    max(abs(err))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% Ivan Selesnick
% selesi@poly.edu

N = length(x);

h0 = af{1};
h1 = af{2};
h2 = af{3};

L0 = length(h0)/2;
L1 = length(h1)/2;
L2 = length(h2)/2;

% symmetric extension
A = L0;
xe = [x(A:-1:1) x x(N:-1:N-A+1)];
% lowpass filter
lo = conv(xe, h0);
% extract valid part
lo = lo(2*L0-1+2*[1:N/2]);

% symmetric extension
A = L1;
xe = [x(A:-1:1) x x(N:-1:N-A+1)];
% highpass filter 1
hi1 = conv(xe, h1);
% down-sample and extract valid part
hi1 = hi1(2*L1-2+2*[1:N/2+1]);
% normalize
hi1(1) = hi1(1)/sqrt(2);
hi1(N/2+1) = hi1(N/2+1)/sqrt(2);

% symmetric extension
A = L2;
xe = [x(A:-1:1) x x(N:-1:N-A+1)];
% highpass filter 2
hi2 = conv(xe, h2);
% down-sample and extract valid part
hi2 = hi2(2*L1-2+2*[2:N/2]);
