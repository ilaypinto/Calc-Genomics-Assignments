%%%%%%%%%%%%%%%%%%%%%%%%%
%% HW5 - Calc Genomics %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the main script for HW 5 in Calc Genomics course, 
% Questions no 1,4,5
% Make sure 'Yeast_ORFS.mat' is in the same file as script

clear all; clc; close all;

%% Q1 - Calculating starting probabilities and transition matrix for sequence
disp('Q1')
disp('__________')

% Define needed parameters
sequence = 'GTAAAGCACAAAAAAAAAAGGGGAATTATTTTGGGGACGGGCAAACCCAAGGGAACAAAATTTAAAAAAAATTTCTAAACTCCAG';
L = length(sequence);

% Count the number of Nucletides
A_count = sum(sequence=='A');C_count = sum(sequence=='C'); 
T_count = sum(sequence=='T');G_count = sum(sequence=='G');

% Computing starting probabilities
A_start = A_count/L; C_start = C_count/L;
T_start = T_count/L; G_start = G_count/L;

starting_probabilities = table(A_start,C_start,G_start,T_start);
disp('Starting Probabilities:')
disp(starting_probabilities)

% Calculate transition probability from A
A_A = length(strfind(sequence,'AA'))/A_count; A_C = length(strfind(sequence,'AC'))/A_count;
A_T = length(strfind(sequence,'AT'))/A_count; A_G = length(strfind(sequence,'AG'))/A_count;

% Calculate transition probability from C
C_A = length(strfind(sequence,'CA'))/C_count; C_C = length(strfind(sequence,'CC'))/C_count;
C_T = length(strfind(sequence,'CT'))/C_count; C_G = length(strfind(sequence,'CG'))/C_count;

% Calculate transition probability from T
T_A = length(strfind(sequence,'TA'))/T_count; T_C = length(strfind(sequence,'TC'))/T_count;
T_T = length(strfind(sequence,'TT'))/T_count; T_G = length(strfind(sequence,'TG'))/T_count;

% Calculate transition probability from G
G_A = length(strfind(sequence,'GA'))/G_count; G_C = length(strfind(sequence,'GC'))/G_count;
G_T = length(strfind(sequence,'GT'))/G_count; G_G = length(strfind(sequence,'GG'))/G_count;

% Define the transition matrix
A = [A_A;C_A;G_A;T_A]; % Define A row
C = [A_C;C_C;G_C;T_C]; % Define C row
T = [A_T;C_T;G_T;T_T]; % Define T row
G = [A_G;C_G;G_G;T_G]; % Define G row

From_To = ['A';'C';'G';'T']; % Row axis
transition_matrix = table(From_To,A,C,G,T);
disp('Transition Matrix:')
disp(transition_matrix)

%% Q4,5 - Viterbi Algorithm
disp('Q4&5')
disp('_________')
% Q4 - myviterbi
% See 'myviterbi' function under 'Local Functions' section in the end of
% this script
disp('Q4 - See "myviterbi" function under Local Functions')

% Q5 - Most likely hidden sequence for ORF no.58 as observation

% Define relevant parameters
load('Yeast_ORFS.mat')        % Load ORF data
Nucleotides = 'ACGT';         % Define Nucleotides
Hidden_States_name = 'HL';    % Define Hidden States
Hidden_States_num = [1,2];
initPr = [0.5, 0.5];          % Define starting probabilities for hidden states
transitionPr = [0.5, 0.5; 0.4, 0.6]; % Define transition matrix
emissionPr = [0.2, 0.3, 0.3, 0.2; 0.3, 0.2, 0.2, 0.3]; % Define emission matrix
observations = zeros(1,length(orfs{58})); % Observations is ORF no.58

% Find the position of each nucleotide in observed ORF
for i = 1:length(Nucleotides) 
    observations(strfind(orfs{58}, Nucleotides(i))) = i;
end

% Find the most likely scenario using 'myviterbi' in Local Functions
most_likely_states_num = myviterbi(Hidden_States_num,initPr, transitionPr, emissionPr, observations);
most_likely_states_num = int2str(most_likely_states_num);

% Create an 'H' & 'L' hidden state vector according to 
% the most likely scenario calculated earlier
for i = 1:length(Hidden_States_name)
    most_likely_states(strfind(most_likely_states_num, int2str(i))) = Hidden_States_name(i);
end

most_likely_states_num = str2num(most_likely_states_num); % Turn back to numbers

disp('Q5:')
disp('For the Hidden States in H&L form, see "most_likely_states"')
disp('For the Hidden States in numeric form, see "most_likely_states_num"')

%% Local Functions

function most_likely = myviterbi(hidden_states,initPr,transitionPr,emissionPr,observations)

% Define Parameters for calculations
% Probabilities are in log as they are very small
    most_likely = zeros(1, length(observations));
    probabilities = zeros(length(hidden_states), length(observations));
    probabilities(:,1) = log(initPr)' + log(emissionPr(:,observations(1)));
    pointer = zeros(length(hidden_states), length(observations));
    
    for i = 2:length(observations)
        for j = 1:length(hidden_states)
            % Add the transition to the mix
            temp_probability = probabilities(:,i-1) + log(transitionPr(:,j));
            [posibilities, posibility_location] = max(temp_probability);
            max_probability = posibilities(1); pointer(j,i) = posibility_location(1);
            % Add the emission to the mix
            probabilities(j,i) = max_probability + log(emissionPr(j,observations(i)));
        end
    end
    
% Most likely scenario - max probability!
    most_likely(end) = find(probabilities(:,end) == max(probabilities(:,end)));
    
% Backtrack using the pointer
    for i = 1:length(most_likely)-1
        most_likely(end-i) = pointer(most_likely(end-i+1),end-i+1);
    end
end