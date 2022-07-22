%% Q2 Genomics
% A function, solving a gene(as defined in this section)
% string into its sequence
function aa_seq=nt_2_aa(nt_vec)
    % computing only if the sequence give in in pairs
    % otherwise, we assume an error occurred in data aquisition

    if mod(length(nt_vec),2)==0
        L=floor(length(nt_vec)/2);

        % define a string to build on, as apending as is gives ASCII values
        % of the characters
        aa_seq='S';

        % solver
        for i=1:L
            codon=[nt_vec(i*2-1) nt_vec(i*2)];
            if strcmp(codon, 'aa')
                aa_seq=strcat(aa_seq,'X');
            elseif strcmp(codon, 'ab')
                aa_seq=strcat(aa_seq,'Y');
            elseif strcmp(codon, 'ba')
                aa_seq=strcat(aa_seq,'Z');
            elseif strcmp(codon,'bb')
                aa_seq=strcat(aa_seq,'W');
            end

            % removing the starter string
            aa_seq=strip(aa_seq,'S');
        end
    else
        % error message
        disp('Given gene length is odd, therefore could not be translated')
        
    end
end