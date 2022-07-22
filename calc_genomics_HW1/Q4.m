%% Q4 Genomics
clear all; clc; close all;

%% Subsection 1-2

% uploading E_coli_ORF.mat and computing number of appearences for
% each gene
data = load('E_cooli_ORF.mat');
data = data.gene_ORF;

% All possible codons
codons=["ATT" "ATC" "ATA" "CTT" "CTC" "CTA" "CTG" "TTA" "TTG" "GTT" "GTC"...
    "GTA" "GTG" "TTT" "TTC" "ATG" "TGT" "TGC" "GCT" "GCC" "GCA" "GCG" "GGT"...
    "GGC" "GGA" "GGG" "CCT" "CCC" "CCA" "CCG" "ACT" "ACC" "ACA" "ACG" "TCT"...
    "TCC" "TCA" "TCG" "AGT" "AGC" "TAT" "TAC" "TGG" "CAA" "CAG" "AAT" "AAC"...
    "CAT" "CAC" "GAA" "GAG" "GAT" "GAC" "AAA" "AAG" "CGT" "CGC" "CGA" "CGG"...
    "AGA" "AGG" "TAA" "TAG" "TGA"];

% finding the number of codons inside each gene
codons_in_genes=zeros(length(data),length(codons));
for i=1:length(data)
   for j=1:floor(length(data{i})/3)
        instance=find(data{i}(j*3-2:j*3)==codons);
        codons_in_genes(i,instance)=codons_in_genes(i,instance)+1;
   end
end
%% Subsection 3

% finding the most abundant codon in each gene, and naming him, arranged
% in a cell array
[~,I]=max(codons_in_genes,[],2);
abundant_codons={};
for i=1:length(I)
    abundant_codons{i}=codons(I(i));
end
abundant_codons = abundant_codons';
%% Subsection 4

% finding the longest gene and creating a bar chart, numbering all the
% codons in it
length_vec=zeros(length(data),1);
for i=1:length(data)
    length_vec(i)=length(data{i});
end
[~,I]=max(length_vec); yaxis=codons_in_genes(I,:);
xaxis = categorical({'ATT' 'ATC' 'ATA' 'CTT' 'CTC' 'CTA' 'CTG' 'TTA' 'TTG' 'GTT' 'GTC'...
    'GTA' 'GTG' 'TTT' 'TTC' 'ATG' 'TGT' 'TGC' 'GCT' 'GCC' 'GCA' 'GCG' 'GGT'...
    'GGC' 'GGA' 'GGG' 'CCT' 'CCC' 'CCA' 'CCG' 'ACT' 'ACC' 'ACA' 'ACG' 'TCT'...
    'TCC' 'TCA' 'TCG' 'AGT' 'AGC' 'TAT' 'TAC' 'TGG' 'CAA' 'CAG' 'AAT' 'AAC'...
    'CAT' 'CAC' 'GAA' 'GAG' 'GAT' 'GAC' 'AAA' 'AAG' 'CGT' 'CGC' 'CGA' 'CGG'...
    'AGA' 'AGG' 'TAA' 'TAG' 'TGA'});
xaxis = reordercats(xaxis,{'ATT' 'ATC' 'ATA' 'CTT' 'CTC' 'CTA' 'CTG' 'TTA' 'TTG' 'GTT' 'GTC'...
    'GTA' 'GTG' 'TTT' 'TTC' 'ATG' 'TGT' 'TGC' 'GCT' 'GCC' 'GCA' 'GCG' 'GGT'...
    'GGC' 'GGA' 'GGG' 'CCT' 'CCC' 'CCA' 'CCG' 'ACT' 'ACC' 'ACA' 'ACG' 'TCT'...
    'TCC' 'TCA' 'TCG' 'AGT' 'AGC' 'TAT' 'TAC' 'TGG' 'CAA' 'CAG' 'AAT' 'AAC'...
    'CAT' 'CAC' 'GAA' 'GAG' 'GAT' 'GAC' 'AAA' 'AAG' 'CGT' 'CGC' 'CGA' 'CGG'...
    'AGA' 'AGG' 'TAA' 'TAG' 'TGA'});
figure;
bar(xaxis,yaxis); title(['Longest Gene Codon Distribution (gene No.', num2str(I),')']);
xlabel('Codons'); ylabel('No. of appearances in Gene');

