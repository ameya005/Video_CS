function w = ddwt(x, J, af)

% Double-Density Wavelet Transform
%
% USAGE:
%    w = ddwt(x, J, af)
% INPUT:
%    x - N-point vector with N divisible by 2^J
%    J - number of stages
%    af - analysis filters (even length)
%    af{i} - filter i (i = 1,2)
% OUTPUT:
%    w{j}{i} - wavelet coefficients (j = 1..J, i = 1,2)
%    w{J+1} - scaling coefficients
% EXAMPLE:
%    [af, sf] = filters1;
%    x = rand(1,128);
%    w = ddwt(x,3,af);
%    y = ddwti(w,3,sf);
%    err = x - y; 
%    max(abs(err))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% Ivan Selesnick
% selesi@poly.edu

for j = 1:J
    [x w{j}{1} w{j}{2}] = afb(x, af);
end
w{J+1} = x;
