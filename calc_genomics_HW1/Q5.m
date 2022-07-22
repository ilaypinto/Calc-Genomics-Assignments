%% Q5 Genomics
clear all; clc; close all;

%% Subsection 1-3

% uploading 'yeast_genes.fa', and coverting to RNA and Amino Acids
% make sure 'yeast_genes.fa' is in the same file as script
[h,s]=fastaread('yeast_genes.fa.txt');
RNA = cellfun(@dna2rna,s,'UniformOutput',false);
AA = cellfun(@nt2aa,s,'UniformOutput',false);
%% Subsection 4

% Length(in number of codons) of each sequence
L = cellfun(@length,s,'UniformOutput',false);
L = cell2mat(L)/3;
%% Subsection 5-6

% using 'CG_content_calc' function on the sequence.
% make sure that the function is in the same file as this script
GC_content = CG_content_calc(s);
%% Subsection 7

% creating a cell array, rearranging all found information in previous
% sections into a single array
all_info(:,1) = h';
all_info(:,2) = s';
all_info(:,3) = AA';
all_info(:,4) = num2cell(L');
all_info(:,5) = GC_content';
