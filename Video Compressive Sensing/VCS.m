%%%%%%INPUT%%%%%
%Run inputImgs to generate the samples Matrix. This can be loaded everytime
samples = load('samples_first.mat');
addpath('./Framelet2X/ToolBox');    % For ddwt and iddwt - framelt functions
X = samples.X;  %input images- vectorized
tr = samples.tr; %input images s 
y = samples.y;
S = samples.S;
n = samples.n; % num of images;
m = samples.m; % num of samples taken totally from video volume
phi = samples.phi; %Measurement matrix
J = 50;
%input video volume (before compressive measuremeent)
%X = ; %n x J - Surveillance Video = [x_1 x_2 .... x_J] where x_j is (n x 1)

%compressive measurement matrix
%phi = ; %m x n - Permuted Walsh Hadamard matrix - there are J such phi which are randomly generated

%compressive measurement / input
% y = zeros(m,1);
% for i=1:J
%     phi = generate_random_phi();
%     y = y + phi*X(:,i); %m x 1
%     if i==1             %stores J phi's 
%         store = phi;    
%     else
%         store = [store phi];
%     end
% end

%sparse basis for recovery
% for j=1:J
%     W = generate_wavelet_transform();
%     if j==1         %append all wavelet transforms for all J frames
%         W1=W; 
%     else
%         W1 = [W1 W];
%     end
% end
[af, sf] = filters1;  %Filters for ddwt

%%%%%%INITIALIZE%%%%%%
%initialize the following
X1 = zeros(n, J) ; %n x J - stationary background - low rank = [x1_1 x1_2 .... x1_J] where x1_j is (n x 1)
X2 = zeros(n, J); %n x J - moving object - sparse = [x2_1 x2_2 .... x2_J] where x2_j is (n x 1)
Z1 = ones(n, J); %n x J - splitting variable of X1 = [z1_1 z1_2 .... z1_J] where z1_j is (n x 1)
Z2 = ones(3*n/2, J) ; %2n x J - splitting variable of (W1 o X1) = [z2_1 z2_2 .... z2_J] where z2_j is (2n x 1)
Z3 = ones(3*n/2, J); %2n x J - splitting variable of (W2 o X2) = [z3_1 z3_2 .... z3_J] where z3_j is (2n x 1)
lambda1 = zeros(n, J); %n x J - lagrangian multiplier for (X1 - Z1) = [lambda1_1 ... lambda1_J] where lambda1_j is (n x 1)
lambda2 = zeros(3*n/2, J); %2n x J - lagrangian multiplier for (W1 o X1 - Z2) = [lambda2_1 ... lambda2_J] where lambda2_j is (2n x 1)
lambda3 = zeros(3*n/2, J); %2n x J - lagrangian multiplier for (W2 o X2 - Z3) = [lambda3_1 ... lambda3_J] where lambda3_j is (2n x 1)
lambda4 = zeros(m,1); %m x 1 - lagrangian multiplier for [phi * (X1+X2) - y] = [lambda4_1 ... lambda4_m] where lambda4_k is (1 x 1)
beta1 = 0.25; %1 x 1 - penalty parameter for (X1 - Z1)
beta2 = 0.25; %1 x 1 - penalty parameter for (W1 o X1 - Z2)
beta3 = 0.25; %1 x 1 - penalty parameter for (W2 o X2 - Z3)
beta4 = 0.25; %1 x 1 - penalty parameter for [phi * (X1+X2) - y]
gamma = 0.01; %1 x 1 - step size for algorithm convergence
mu1 = 0.01; %1 x 1 - minimization coefficient of Z1
mu2 = 0.01; %1 x 1 - minimization coefficient of Z2
mu3 = 0.01; %1 x 1 - minimization coefficient of Z3
%% tau = ; %1 x 1 - threshold parameter for SVT 
converged = 1;

iter = 0;
%%%%%%MINIMIZATION-1%%%%%%
%%%compute X1,X2 by holding Z1,Z2 and Z3 as constant, from nabla L_A%%%
while iter~=10
    iter = iter+1
    M1 = zeros(n,J);
    M2 = zeros(n,J);
%%%calculate M1 matrix%%%
%     for i=1:n
%         i
%         for j=1:J
%             phi_j = phi(:,:,j);  %select required phi from stored phi's
%             for k=1:m
%                 M1(i,j) = M1(i,j) + lambda4(k)*phi_j(k,i);
%                 
%             end
%         end
%     end
    for k=1:m
        M1 = lambda4(k)*reshape(phi(k,:,:), n, J);
    end
    fprintf('M1 done!');
%     %%%calculate M2 matrix%%%
%     for i=1:n
%         for j=1:J
%             phi_j = phi(:,:,j);  %select required phi from stored phi's
%             for k=1:m
%                 M2(i,j) = M2(i,j) + phi_j(k,i)*y(k,1);
%             end
%         end
%     end
    for k=1:m
        M2 = y(k,1)*reshape(phi(k,:,:),n,J);
    end
    M2 = beta4*M2;
    fprintf('M2 done!');

    %%%%%%%MINIMIZATION-1-1%%%%%%%%%
    %%%equating gradient X1 to zero%%%
    % lambda1 + (W1^T o lambda2) + M1 + beta1*(X1-Z1) + beta2*(X1 - W1^T o Z2)
    % +beta4*(sum of columns of X1) - M2 = 0

    %%%or%%%

    % (beta1+beta2)*X1 + beta4*(sum of columns of X1)
    % = -(lambda1 + (W1^T o lambda2) + M1) + beta1*Z1 + beta2*(W1^T o Z1) + M2
    % = constant = C

    alpha = beta1+beta2;
    beta = beta4;

    %%%calculate constant C%%%
    %W1^T = [W1_1^T W1_2^T ..... W1_J^T]
    %C2 = (W1^T o Z1) = [W1_1^T*z1_1 W1_2^T*z1_2 ..... W1_J^T*z1_J]
    C2 = zeros(n,J);
    %C1 = W1^T;
    for j=1:J      %n x J
        C2(:, j) = invFrameTr(Z2(:,j));% C1(:n*(j-1)+1:n*j)*Z1(:,j); % n x 1 %%%%%perform Inverse DDWT over here <-- insert function
    end
    C2 = beta2*C2;
    fprintf('C2 done!');
    %C3 = W1^T o lambda2
    C3 = zeros(n,J);
    for j=1:J      %n x J
        C3(:,j) = invFrameTr(lambda2(:,j));%C1(:n*(j-1)+1:n*j)*lambda2(:,j); % n x 1  %%%%%perform Inverse DDWT over here <-- insert function
    end
    fprintf('C3 done!');
    C = -lambda1 -C3 - M1 + beta1*Z1 + C2 + M2; % n x J

    %solution to alpha*X1 + beta*(sum of columns of X1) = C can be obtained by
    %using following inversion formula
    %Ax_k = c_k --> x_k = (A^-1)*c_k
    %A = [alpha+beta beta beta.....;beta alpha+beta beta.....;......;....;beta beta....alpha+beta]
    alpha = alpha+beta;
    %A = [alpha beta beta.....;beta alpha beta.....;......;....;beta beta....alpha]
    A = zeros(J,J);
    A = alpha*ones(J,J) + beta*eye(J);
    fprintf('A done!');
    for k=1:n
        X1(k,:) = C(k,:)*(A^-1)'; %final X1 from minimization
    end
    fprintf('X1 done!');
    %%%%%%%%MINIMIZATION-1-2%%%%%%%%%%%
    %%%equating gradient X1 to zero%%%
    % (W2^T o lambda3) + M1 + beta3*(X2 - W2^T o Z3)
    % + beta4*(sum of columns of X1) - M2 = 0

    %%%or%%%

    % (beta3)*X2 + beta4*(sum of columns of X1)
    % = -((W2^T o lambda3) + M1) + beta3*(W2^T o Z3) + M2
    % = constant = C

    a = beta3;
    b = beta4;

    %%%calculate constant C%%%
    %W2^T = [W2_1^T W2_2^T ..... W2_J^T]
    %D2 = (W2^T o Z1) = [W2_1^T*z2_1 W2_2^T*z2_2 ..... W2_J^T*z2_J]
    D2 = zeros(n,J);
    %D1 = W2^T;
    for j=1:J      %n x J
        D2(:,j) = invFrameTr(Z2(:,j));%D1(:n*(j-1)+1:n*j)*Z2(:,j); % n x 1  %%%%%perform Inverse DDWT over here <-- insert function
    end
    D2 = beta3*D2;
    %C3 = W1^T o lambda2
    D3 = zeros(n,J);
    for j=1:J      %n x J
        D2(:,j) = invFrameTr(lambda3(:,j));%D1(:n*(j-1)+1:n*j)*lambda3(:,j); % n x 1  %%%%%perform Inverse DDWT over here <-- insert function
    end
    D = -D3 - M1 + D2 + M2; % n x J
    fprintf('D done!');
    %solution to a*X2 + b*(sum of columns of X2) = D can be obtained by
    %using following inversion formula
    %Mw_k = d_k --> w_k = (M^-1)*d_k
    %M = [alpha beta beta.....;beta alpha beta.....;......;....;beta beta....alpha]
    M = zeros(J,J);
    M = alpha*ones(J,J) +(beta-alpha)*eye(J);   
%     for i=1:J
%         for j=1:J
%             if j==i
%                 M(i,j) = alpha;
%             else
%                 M(i,j) = beta;
%             end
%         end
%     end
    
    fprintf('M done!');
    for k=1:n
        %X2(k,:) = (M^-1)*D(k,:); %final X2 from minimization
        X2(k,:) = D(k,:) * (M^-1)';
    end
    fprintf('X2 done!');
    %%%%%%MINIMIZATION-2%%%%%
    %%%compute Z1,Z2,Z3 by Singular Value Thresholding and Shrinkage Formula%%%

    %%%Z1 computation%%%%
    %%%SVT function%%%%
    tau = mu1/beta1;
    X = X1 - lambda1/beta1;
    Z1 = SVT(tau,X,J);

    %%%Z2,Z3 computation%%%%
    %%%Shrinkage operator%%%
    tau2 = mu2/beta2;
    tau3 = mu3/beta3;
    x2 = zeros(3*n/2,J);
    x3 = zeros(3*n/2,J);
    for i=1:J
        x2(:,i) = frameTr(X1(:,i));%%%% apply DDWT here <-- write the function
    end
    for i=1:J
        x2(:,3) = frameTr(X2(:,i));%%%% apply DDWT here <-- write the function
    end
    Z2 = Shrinkage(x2- (lambda2/beta2),tau2,J);
    Z3 = Shrinkage(x3- (lambda3/beta3),tau3,J);
    fprintf('X1,X2,Z1,Z2,Z3 done!');
    %Updatng lambda
    lambda1 = lambda1 - (gamma*beta1)*(X1-Z1);
    lambda2 = lambda2 - (gamma*beta2)*(x2 - Z2);
    lambda3 = lambda3 - (gamma*beta3)*(x3 - Z3);
    phi_x = zeros(m,1);
    X_full = X1+X2;
    for j=1:J
        phi_j = phi(:,:,j);
        phi_x = phi_x + phi_j*X_full(:,j);
    end    
    lambda4 = lambda4 - (gamma*beta4)*(phi_x - y);
    fprintf('Lambdas done!');
    fprintf('Rank: %d', rank(X1));
    fprintf('norm: %d', norm(X1,1));
    %convergence check
    if(iter >= 10)
        converged = 0;
    end
end                             %end of while