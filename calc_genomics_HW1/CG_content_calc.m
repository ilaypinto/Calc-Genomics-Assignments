function [GC_content] = CG_content_calc(sequence)
% Computing GC[%] as stated by the formula given, per cell given in a cell
% array

% acquire N values using 'basecount' 
N=cellfun(@basecount,sequence,'UniformOutput',false);

% compute GC[%] according to formula(see local functions)
GC_content = cellfun(@GC_mid_calc,N,'UniformOutput',false);
end
%% Local Functions
function [GC_content] = GC_mid_calc(N)
GC_content = ((N.G + N.C)/(N.G + N.C + N.A + N.T))*100;
end