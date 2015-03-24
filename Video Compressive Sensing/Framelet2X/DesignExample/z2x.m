function p = z2x(h)
% p = z2x(h)
% Implements the change of variables
% x = (-z + 2 - 1/z)/4
% where h(z) is a odd-length symmetric filter

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic University
% Brooklyn, NY, USA

N = length(h);
M = (N-1)/2;
p = [];
g = 1;
for k = 1:M
   g = conv(g,[-1 2 -1]/4);
end
for k = 0:M
   [q,r] = deconv(h,g);
   p(M+1-k) = q;
   h = r(2:end-1);
   g = deconv(g,[-1 2 -1]/4);
end
p = p(end:-1:1);
