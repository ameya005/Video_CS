function [X, Tr, S, y, phi, n, m] = inputImgs(path_file)

    cd(path_file)
    pwd
    files = dir('.');
    img = imread(files(3).name);
    img = imresize(img, [64, 128]);
    [s, t] = size(rgb2gray(img));
    n = s*t;
    X = zeros(n,1);
    Tr = zeros(n,50); %Walsh Hadamard creates a matrix of the closest 2's power
    m = ceil(0.04*n); % percentage of samples ---vary this to get different results
    S = zeros(m,50);
    phi = zeros(m,n, 50);
    tempI = eye(s*t);
    %phi_j = reshape(tempI,n,1);
    temptr = fwht(tempI);
    for i = 3:52
        img = imread(files(i).name);
        img = rgb2gray(img);
        img = imresize(img, [64, 128]);
        x = reshape(img, n, 1);
        tr = fwht(double(x));
        %size(tr)
        X(:,i-2) = x;
        Tr(:,i-2) = tr;
    end
    phi_j = zeros(m, n);
    %size(X)
    r = randi(n, 1, m);
    for j=1:50
        for i=1:m
            S(i,j) = Tr(r(i),j);
            phi_j(i,:) = temptr(r(i), :);
        end
        phi(:,:,j) = phi_j(:,:);
    end
    
    y = sum(S, 2);
    cd ..
    save('samples_first.mat');  % creates a sample matrix;
    %X(96,:)
    %imshow(X,[]);
    
end
        