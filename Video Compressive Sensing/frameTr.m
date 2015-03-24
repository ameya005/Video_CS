function W = frameTr(X)

    addpath('./Framelet2X/ToolBox');
    [af, sf] = filters1;
    [s1, s2] = size(X);
    cellW = ddwt(X', 1, af);
    W = [ cellW{1}{1}, cellW{1}{2}, cellW{2}];
    W = W';
end
    