function [ D ] = SVT( tau,X,J )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[U S V] = svd(X,0);
th = S-tau*eye(J,J);
for j=1:J
    if th(j,j) < 0
        th(j,j) = 0;
    end
end

D = U*th*V';
end

