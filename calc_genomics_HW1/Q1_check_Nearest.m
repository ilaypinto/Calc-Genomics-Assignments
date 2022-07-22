%% Check Code fore Q1
% checks it works for both matrices and vectors
findNearest(rand(5),1)
findNearest(rand(5,1),1)

% check it gives all values, when matrix is uniform
findNearest(ones(5),1)

% check it gives the right ones only for 2 indices
k=ones(5);k(1,1)=0;k(3,4)=0;
findNearest(k,0)
findNearest(k,1)