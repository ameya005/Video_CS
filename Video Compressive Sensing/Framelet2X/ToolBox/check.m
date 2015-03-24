%% Verifying MATLAB Programs for the Symmetric Double-Density DWT

%% Check perfect reconstruction property of afb and sfb

disp('Check PR property of afb and sfb:')
[af, sf] = filters1;
x = rand(1,64);
[lo, hi1, hi2] = afb(x, af);
y = sfb(lo, hi1, hi2, sf);
err = x - y;
MaxErr = max(abs(err));
fprintf('Maximum reconstruction error = %f\n\n',MaxErr);

%% Check Parseval's energy condition

disp('Check Parseval''s energy condition:')
E = x*x';
E2 = lo*lo' + hi1*hi1' + hi2*hi2';
fprintf('Energy of signal = %f\n', E);
fprintf('Energy of subband signals = %f\n', E2);
fprintf('The difference = %f\n\n', E-E2);

%% Check perfect reconstruction property of ddwt and ddwti

disp('Check PR property of ddwt and ddwti:')
[af, sf] = filters1;
x = rand(1,128);
w = ddwt(x,3,af);
y = ddwti(w,3,sf);
err = x - y;
MaxErr = max(abs(err));
fprintf('Maximum reconstruction error = %f\n\n',MaxErr);

%% Display filters: impulse responses and frequency responses

h0 = af{1};
h1 = af{2};
h2 = af{3};
N0 = length(h0);
N1 = length(h1);
N2 = length(h2);

% compute frequency responses:
L = 2^10;
H0 = fft(h0,L);
H1 = fft(h1,L);
H2 = fft(h2,L);
H0 = H0(1:L/2+1);
H1 = H1(1:L/2+1);
H2 = H2(1:L/2+1);
f = [0:L/2]*2/L;

figure(1)
clf
subplot(4,2,1)
stem(0:N0-1,h0,'.')
title('h0(n)')
ax = [-1 N1 -1 1];
axis(ax);

subplot(4,2,2)
plot(f,abs(H0));
title('H0(\omega)')

subplot(4,2,3)
stem(0:N1-1,h1,'.')
title('h1(n)')
axis(ax);

subplot(4,2,4)
plot(f,abs(H1));
title('H1(\omega)')

subplot(4,2,5)
stem(0:N2-1,h2,'.')
title('h2(n)')
xlabel('n')
axis(ax);

subplot(4,2,6)
plot(f,abs(H2));
title('H2(\omega)')
xlabel('\omega/\pi')

orient tall
print -dpsc check
pause(1)

%% Analyze a test signal with analysis filter bank:

disp('Analyze a test signal with the 3-channel analysis filter bank')
n = 1:58;
x = double(((n > 10)&(n<24))|(n>41));
[af, sf] = filters1;
[lo, hi1, hi2] = afb(x, af);

figure(1)
clf
subplot(4,1,1)
stem(x,'.')
title('TEST SIGNAL')

subplot(4,1,2)
stem(lo,'.')
title('LOW-PASS OUTPUT')

subplot(4,1,3)
stem(hi1,'.')
title('BAND-PASS OUTPUT')

subplot(4,1,4)
stem(hi2,'.')
title('HIGH-PASS OUTPUT')

orient tall
print -append -dpsc check
pause(1)

%% Analyze a test signal with double-density wavelet transform:

disp('Analyze a test signal with the double-density DWT')
x = [zeros(1,19) [1:48]/48 zeros(1,21) ones(1,40)];
N = length(x);
[af, sf] = filters1;
w = ddwt(x, 3, af);

figure(1)
clf
subplot(7,1,1)
stem(x,'.')
title('TEST SIGNAL')
axis tight off

subplot(7,1,2)
stem(w{1}{1},'.')
title('w\{1\}\{1\}')
axis tight off

subplot(7,1,3)
stem(w{1}{2},'.')
title('w\{1\}\{2\}')
axis tight off

subplot(7,1,4)
stem(w{2}{1},'.')
title('w\{2\}\{1\}')
axis tight off

subplot(7,1,5)
stem(w{2}{2},'.')
title('w\{2\}\{2\}')
axis tight off

subplot(7,1,6)
stem(w{3}{1},'.')
title('w\{3\}\{1\}')
axis tight off

subplot(7,1,7)
stem(w{3}{2},'.')
title('w\{3\}\{2\}')
axis tight off

orient tall
print -append -dpsc check
pause(1)

%% Plot the first wavelet

N = 512;
n = 0:N-1;
z = zeros(1,N);
[af, sf] = filters1;
J = 4;
wz = ddwt(z, J, af);
w = wz;
L = N/2^J;
% compute a wavelet that does not meet the signal boundary:
w{J}{1}(L/2) = 1;
y1 = ddwti(w, J, sf);
w = wz;
% compute wavelets that do meet the signal boundary:
w{J}{1}(L-1) = 1;
y2 = ddwti(w, J, sf);
w = wz;
w{J}{1}(L) = 1;
y3 = ddwti(w, J, sf);
w = wz;
w{J}{1}(L+1) = 1;
y4 = ddwti(w, J, sf);

figure(1)
clf
subplot(4,1,1)
plot(n,y1)
title('FIRST WAVELET')
xtight

subplot(4,1,2)
plot(n,y2)
xtight

subplot(4,1,3)
plot(n,y3)
xtight

subplot(4,1,4)
plot(n,y4)
xtight

orient tall
print -append -dpsc check
pause(1)

%% Plot the second wavelet

w = wz;
w{J}{2}(L/2) = 1;
y1 = ddwti(w, J, sf);
w = wz;
w{J}{2}(L-3) = 1;
y2 = ddwti(w, J, sf);
w = wz;
w{J}{2}(L-2) = 1;
y3 = ddwti(w, J, sf);
w = wz;
w{J}{2}(L-1) = 1;
y4 = ddwti(w, J, sf);

figure(1)
clf
subplot(4,1,1)
plot(n,y1)
title('SECOND WAVELET')
xtight

subplot(4,1,2)
plot(n,y2)
xtight

subplot(4,1,3)
plot(n,y3)
xtight

subplot(4,1,4)
plot(n,y4)
xtight

orient tall
print -append -dpsc check
pause(1)
