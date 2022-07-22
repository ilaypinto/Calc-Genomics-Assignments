%%%%%%%%%%%%%%%%%%%%%%%
%% Calc Genomics HW 6%%
%%%%%%%%%%%%%%%%%%%%%%%

% This is the main script for calc genomics HW 6.
% Make sure "yeast_genes.fa" and "MTDR_meiosis_yeast"
% are in the same file as this script before running

clear all; close all; clc;

%% Q1 - GC content, Hierarchial Clustering
disp('Q1 - GC content, Hierarchial Clustering')
disp('_______________________________________')
% Section 1
[~,s] = fastaread('yeast_genes.fa.txt'); % Load Data
GC = CG_content_calc(s(:,1:5));          % Calc GC content for first 5
GC = cell2mat(GC);
disp('Section 1 - Initial GC content calculated and is')
disp(GC)

% Section 2 - Euclidean Distance
distance = [];
for i = 1:5
     distance = [distance,Distance(GC(i),GC)]; % Calc distance
end
disp('Section 2 - Euclidean Distance calculated and is')
disp(distance)

% Section 3-4 - Hierarchial Clustering
disp('Section 3 - see manual calculation in answers pdf')
link = linkage(distance,'complete');
figure;
dendrogram(link); title('GC content Hierarchial Clustering');
disp('Section 4 - Dendrogram computed!')

%% Q2 - Cluster Manual Computation

disp('Q2 - Cluster Manual Computation')
disp('_______________________________')
disp('See manual calculation in answers pdf')

%% Q3 - MTDR, Zscore, etc.
disp('Q3 - MTDR, Zscore, etc.')
disp('_______________________')

load('MTDR_meiosis_yeast.mat') % Load MTDR

% Section 1 - Z score
Z_score = ZZZ(ER_matrix); % calc Z score
disp('Section 1 - Z score computed! See Z_score in workspace')

% Section 2 - kmeans handle
kmeans_handle = @(mat,k) kmeans(mat,k,'Distance','correlation','MaxIter',1000);
disp('Section 2 - kmeans handle computed! See kmeans_handle in workspace')

% Section 3 - find optimal k value for Z score
DB_eval = evalclusters(Z_score,kmeans_handle,'DaviesBouldin','KList',1:10);
k_DB = DB_eval.OptimalK;                % Optimal Davies-Bouldin k
Silouette_eval = evalclusters(Z_score,kmeans_handle,'silhouette','KList',1:10);
k_Silouette = Silouette_eval.OptimalK;  % Optimal Silouette k
disp('Section 3 - Optimal k value for Z score were computed:')
disp(['By using Davies-Bouldin, optimal k is ', num2str(k_DB)]);
disp(['By using silhouette, optimal k is ', num2str(k_Silouette)]);

% Section 4 - Plot the clusters

% DB
figure;
index = kmeans_handle(Z_score,k_DB);     
plot(find(index == 1),Z_score(index == 1,1),'r.','MarkerSize',8)
hold on;
plot(find(index == 2),Z_score(index == 2,1),'b.','MarkerSize',8)
plot(find(index == 3),Z_score(index == 3,1),'g.','MarkerSize',8)
hold off;
legend('Cluster 1','Cluster 2','Cluster 3');title('Davies Bouldin');

% Silhouette
figure;
index = kmeans_handle(Z_score,k_Silouette);     
plot(find(index==1),Z_score(index==1,1),'r.','MarkerSize',8)
hold on
plot(find(index==2),Z_score(index==2,1),'b.','MarkerSize',8)
hold off;
legend('Cluster 1','Cluster 2'); title('Silhouette');
disp('Section 4 - See relevant plots')

%% Q4 - Good measurement to evaluate clustering

disp('Q4 - Good measurement to evaluate clustering')
disp('____________________________________________')
disp('See answer for this question in answers pdf')

%% Local Functions

function Z_score = ZZZ(MTDR)

Z_score = zeros(size(MTDR));

    for i = 1:size(MTDR,2)
        for j = 1:size(MTDR,1)
            no_j = MTDR([1:j-1,j+1:end],i);
            Z_score(j,i) = (MTDR(j,i) - mean(no_j))/std(no_j);
        end
    end
end
% GC_content
function GC_content = CG_content_calc(sequence)

N=cellfun(@basecount,sequence,'UniformOutput',false);
GC_content = cellfun(@GC_mid_calc,N,'UniformOutput',false);
end
function [GC_content] = GC_mid_calc(N)
GC_content = ((N.G + N.C)/(N.G + N.C + N.A + N.T))*100;
end

% Euclidean Distance 
function distance = Distance(GCj,GCi)
distance = zeros(length(GCi),1);
    for i=1:length(GCi)
        distance(i) = abs(GCi(i)-GCj);
    end
end
