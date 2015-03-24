function x = invFrameTr(w)

    addpath('./framelet2X/toolBox/');
    [s1 s2] = size(w);
    idx1=s1/3;
    [af, sf] = filters1;
    w_tr{1}{1} = w(1:idx1+1,1)';
  
    w_tr{1}{2} = w(idx1+2:2*idx1, 1)';
    w_tr{2} = w(((2*idx1)+1): s1,1)';
   
    x = ddwti(w_tr, 1, sf);
end
    
