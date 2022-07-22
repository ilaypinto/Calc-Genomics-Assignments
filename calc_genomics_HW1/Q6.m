%% Q6

% finding longest polyA in a HIV sequence from NCBI
% accession number: NC_001802.1
clear all; clc; close all;

%% Q6c-d

% uploading the data and using 'FindPolyA' function.
% make sure the function and data is in the same file as this script
S = genbankread('sequence.gb');
HIV_sequence = S.Sequence;
indices = FindPolyA(HIV_sequence);