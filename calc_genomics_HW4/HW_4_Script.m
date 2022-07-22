%%%%%%%%%%%%%%%%%%%%%%%%
%% Calc Genomics HW 4 %%
%%%%%%%%%%%%%%%%%%%%%%%%

% This is the main script for calc genomics HW4 - Q5
% 'PSSM' Function is in the end of the script under "LOCAL FUNCTIONS"

clear all; clc; close all;

%% Q5a - PSSM Function

% See 'PSSM' under "Local Functions" in the end of this script

%% Q5b - PSSM matrix for given data

load('Yeast_ORFS.mat'); load('Yeast_UTR5.mat'); % Load yeast data
length_UTR5 = 6; length_ORF = 3;                % Define wanted sequence length
Nucleotides = 'ACGT';                           % Nucleotides
start_codon = zeros(4,3);                       % Start codon probabilities check

[A_UTR5,C_UTR5,G_UTR5,T_UTR5] = PSSM(length_UTR5,'UTR5',utrs);
mat_UTR5 = [A_UTR5' C_UTR5' G_UTR5' T_UTR5'];   % Compute PSSM for UTR5s

[A_ORF,C_ORF,G_ORF,T_ORF] = PSSM(length_ORF,'ORF',orfs);
mat_ORF = [A_ORF' C_ORF' G_ORF' T_ORF'];        % Compute PSSM for ORFs

for i = 1:size(orfs,1)                          % Make sure ATG(start codons) 
    ATG_sequence = orfs{i,1}(1:3);              % are in place
    for p = 1:3
        start_codon_finder = strfind(Nucleotides, ATG_sequence(p));
        start_codon(start_codon_finder, p) = start_codon(start_codon_finder, p) +1;
    end
end

start_codon = start_codon/size(orfs,1);   % Find start codons
start_codon = start_codon';               % probability by dividing
                                          % in the number of ORFs  
value = [mat_UTR5; start_codon; mat_ORF]; % Concatenate all calculations

figure;                                   % Figure same as in Tirgul 8
placement = categorical({'-6';'-5';'-4';'-3';'-2';'-1';...
    'A';'T';'G';'1';'2';'3'});
placement = reordercats(placement,{'-6';'-5';'-4';'-3';'-2';...
    '-1';'A';'T';'G';'1';'2';'3'});
bar(placement, value,'stacked'); title('Yeast PSSM - All Genes');
legend('A', 'C', 'G', 'T', 'Location', 'eastoutside');
ylabel('Percent[%]'); xlabel('Nucleotide Position(Relative to ATG)');

%% Q5c - Seqlogo display of the data

my_seqlogo = cell(1, size(orfs,1));
for j = 1:size(orfs,1)                    % Organize sequence for Seqlogo
    my_seqlogo{j} = [utrs{j}(end-5:end),orfs{j}(1:6)];
end
seqlogo(my_seqlogo)

%% Q5d - Finding most probable start codon placement for a given sequence

% Define given sequence
Given_sequence='AAACCGCATGGGATGCGTGCCATGTAAGCGATGTCGTGGAGTATGACTGCCGGGATGATGAATTATGGGCGGTGA';

ATG_indices = strfind(Given_sequence,'ATG'); % Find potential start codons
scores = ones(1, size(ATG_indices,2));     % Define scores for start codons 

for p = 1:size(ATG_indices,2)

    % For each potential ATG, define its 6 nucleotides UTR5
    potential_UTR5=Given_sequence(ATG_indices(p)-6:ATG_indices(p)-1);

    % For each Potential ATG, define its 3 nucleotides ORF
    potential_ORF=Given_sequence(ATG_indices(p)+3:ATG_indices(p)+5);
    
    % Find nucleotide positions for UTR5 and multiply
    for q = 1:length(potential_UTR5) 
        codon_finder = strfind(Nucleotides, potential_UTR5(q));
        scores(p) = scores(p) * mat_UTR5(q,codon_finder);
    end
    
    % Find nucleotide positions for ORF and multiply
    for r = 1:length(potential_ORF)
        codon_finder = strfind(Nucleotides, potential_ORF(r));
        scores(p) = scores(p) * mat_ORF(r,codon_finder);
    end
end

% Find the max score, thus highest probability for the real start codon
real_ATG_index = ATG_indices(scores==max(scores));

%% Local Functions

function [A, C, G, T] = PSSM(wanted_length, type, sequence)

% Define Parameters for PSSM Matrix
Nucleotides = 'ACGT';
probability = zeros(4,wanted_length);

% Cut the relevant part according to data type
for i = 1:size(sequence,1)

    if strcmp(type,'UTR5') % If UTR5, give sequence from end in size of wanted length
        seq_finder = sequence{i,1}(length(sequence{i,1})-wanted_length+1:end);

    elseif strcmp(type,'ORF') % If ORF, skip ATG and give sequence in size of wanted length
        seq_finder = sequence{i,1}(4:3+wanted_length);

    else
        warning('Type was not specified as UTR5 or ORF')
    end

    % Count the Nucleotides
    for j = 1:wanted_length
        indices = strfind(Nucleotides, seq_finder(j)); % Find codon in position
        probability(indices, j) = probability(indices, j) +1; % Add to counter
    end
end

probability = probability/size(sequence,1); % Devide by number of sequences 
                                            % to achieve probability from 
                                            % counter
% Separate each nucleotide probability                                            
A = probability(1, :); C = probability(2, :);
G = probability(3, :); T = probability(4, :);
end