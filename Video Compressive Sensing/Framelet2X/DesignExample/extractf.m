function f = extractf(h,p);
%
% f = extractf(h,p)
% find f such that h = conv(f,p)
%
% When such an f exists, this function has better
% numerical accuracy that the "deconv" command

% Ivan Selesnick
% selesi@poly.edu
% Polytechnic University
% Brooklyn, NY, USA

p = p(:);
h = h(:);
Np = length(p);
Nh = length(h);

C = convmtx(p,Nh-Np+1);

f = C\h;

% check accuracy of result:

SN = 0.000001;     % Small Number
e = max(abs(C*f - h));
% disp(e)
if e > SN
	disp('there is a problem in extracf')
	keyboard
end

f = f';
