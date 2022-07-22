%% Q1 Genomics
% A function which computes the nearest value in array, to a desired value
function ind=findNearest(x, desiredVal)
    if or(size(x,1)==1, size(x,2)==1)
        distance=abs(x-desiredVal);
        min_val=min(distance);
        ind= find(distance==min_val);
    else
        distance=abs(x-desiredVal);
        min_val=min(distance,[],'all');
        index= find(distance==min_val);
        [rows,cols]=ind2sub(size(x),index);
        ind=[rows cols];
    end
end