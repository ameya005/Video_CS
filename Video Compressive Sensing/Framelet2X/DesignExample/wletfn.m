function [w,s,t] = wletfn(h0,h1,K);

% [w,s,t] = wletfn(h0,h1,K);
%
% Computes the scaling function and wavelet
%
% Ivan Selesnick
% selesi@poly.edu
% Polytechnic University
% Brooklyn, NY, USA

if nargin < 3
	K = 7;
end

N0 = length(h0);
N1 = length(h1);

[s,t] = scalfn(h0,K);

L = length(s);

w = sqrt(2)*conv(up(h1,2^(K-1)),s(1:2:L));

L = (N0-1)/2 + (N1-1)/2;

t = [0:2^K*L]/2^K;

w = w(1:(2^K*L+1));


