function rts = dbleroots(p)

% Find roots of a polynomial with double roots:
% P(z) = Q(z)^2.
% Find roots of Q(z).
% Proceed by taking derivative of P(z) to improve
% the numerical accuracy of the root computation.

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic University
% Brooklyn, NY, USA

N = length(p)-1;
p = p(:)';     % ensure p is a row vector
pdiff = p(1:N) .* (N:-1:1);

rts_p = roots(p);
rts_pdiff = roots(pdiff);

rts = rts_p;
for k = 1:N
    [tmp, i] = min(abs(rts(k)-rts_pdiff));
    rts(k) = rts_pdiff(i);
end

trim = zeros(N/2,1);
for i = 1:N/2
    trim(i) = rts(1);
    rts(1) = [];
    [tmp,k] = min(abs(trim(i)-rts));
    rts(k) = [];
end

rts = trim;
