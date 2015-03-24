function [s,t] = scalfn(h,J)
% [s,t] = scalfn(h,J);
% Scaling function obtained by dyadic expansion
% input
%    h : scaling filter
% output
%    s : samples of the scaling function phi(t)
%        for t = k/2^J, k=0,1,2,...
% % Example:
%    h = [1+sqrt(3) 3+sqrt(3) 3-sqrt(3) 1-sqrt(3)]/(4*sqrt(2));
%    [s,t] = scalfn(h);
%    plot(t,s)

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic University
% Brooklyn, NY, USA

if nargin < 2 
   J = 5;
end

N = length(h);
h = h(:).';             % form a row vector

% check sum rules
n  = 0:N-1;
e0 = sum(h) - sqrt(2);
e1 = sum(((-1).^n).*h);
if abs(e0) > 0.0001
   disp('   need: sum(h(n)) = sqrt(2)')
   return
end
if abs(e1) > 0.0001
   disp('   need: sum((-1)^n h(n)) = 0')
   return
end

% Make convolution matrix
H = toeplitz([h zeros(1,N-1)]',[h(1) zeros(1,N-1)]);
% or: H = convmtx(h(:),N);

% Make P matrix
P = sqrt(2)*H(1:2:2*N-1,:);

% Solve for vector
s = [P-eye(N); ones(1,N)] \ [zeros(N,1); 1];
s = s.';                    % phi at integers
L = N;                      % length of phi vector

% Loop through scales
for k = 0:J-1
  s = sqrt(2)*conv(h,s);
  L = 2*L-1;
  s = s(1:L);
  h = up(h,2);
end

% Time axis
t = (0:L-1)*(N-1)/(L-1);
