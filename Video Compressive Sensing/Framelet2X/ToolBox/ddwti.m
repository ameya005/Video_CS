function y = ddwti(w, J, sf)

% Inverse Double-Density Wavelet Transform
%
% USAGE:
%    y = ddwti(w, J, sf)
% INPUT:
%    w - wavelet coefficients
%    J - number of stages
%    sf - synthesis filters
% OUTPUT:
%    y - output signal
% See also ddwt
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

y = w{J+1};
for j = J:-1:1
   y = sfb(y, w{j}{1}, w{j}{2}, sf);
end
