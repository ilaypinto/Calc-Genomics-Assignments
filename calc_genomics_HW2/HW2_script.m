%% Script for relevant questions HW2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
%% Q2: using global and local alignments:%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Q2')
disp('Subsection 1 - Nucleotides')

% Subsection 1- using functions on nucleotides themselves
seq1 = 'CCCGCCGCTGTACGTACGCTAGCTA'; seq2 = 'CCCGCGGCCGTACGCTCACTAGCTCAGTA';
[global_score_NT,global_alignment_NT] = nwalign(seq1,seq2,'Alphabet','NT')
[local_score_NT,local_alignment_NT] = swalign(seq1,seq2,'Alphabet','NT')

disp('Subsection 2 - Amino Acids')

% Subsection 2- translating to amino acids and using the functions again
aa1 = nt2aa(seq1); aa2 = nt2aa(seq2);
[global_score_AA,global_alignment_AA] = nwalign(aa1,aa2,'Alphabet','AA')
[local_score_AA,local_alignment_AA] = swalign(aa1,aa2,'Alphabet','AA')

%% Q5a: Alignments, scores, etc %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Q5a')

%% Subsection 1,2- Global alignment analysis

disp('Subsection 1,2- Global alignment analysis')

% SUBSBACLE Subtilisin Savinase|P29600
sequence1 = 'AQSVPWGISRVQAPAAHNRGLTGSGVKVAVLDTGISTHPDLNIRGGASFVPGEPSTQDGNGHGTHVAGTIAALNNSIGVLGVAPSAELYAVKVLGASGSGSVSSIAQGLEWAGNNGMHVANLSLGSPSPSATLEQAVNSATSRGVLVVAASGNSGAGSISYPARYANAMAVGATDQNNNRASFSQYGAGLDIVAPGVNVQSTYPGSTYASLNGTSMATPHVAGAAALVKQKNPSWSNVQIRNHLKNTATSLGSTNLYGSGLVNAEAATR';

% P41363|ELYA_BACHD Thermostable alkaline protease precursor (Bacillus)
sequence2 = 'MRQSLKVMVLSTVALLFMANPAAASEEKKEYLIVVEPEEVSAQSVEESYDVDVIHEFEEIPVIHAELTKKELKKLKKDPNVKAIEKNAEVTISQTVPWGISFINTQQAHNRGIFGNGARVAVLDTGIASHPDLRIAGGASFISSEPSYHDNNGHGTHVAGTIAALNNSIGVLGVAPSADLYAVKVLDRNGSGSLASVAQGIEWAINNNMHIINMSLGSTSGSSTLELAVNRANNAGILLVGAAGNTGRQGVNYPARYSGVMAVAAVDQNGQRASFSTYGPEIEISAPGVNVNSTYTGNRYVSLSGTSMATPHVAGVAALVKSRYPSYTNNQIRQRINQTATYLGSPSLYGNGLVHAGRATQ';

% Global score and alignment
[g_score_P29600_P41363,g_alignment_P29600_P41363] = nwalign(sequence1,sequence2,'Alphabet','AA','ScoringMatrix','BLOSUM62')

% Alignment Length
g_alignment_P29600_P41363_length = length(g_alignment_P29600_P41363)

% Identity and Similarity Percent
Indentity=sum(g_alignment_P29600_P41363(2,:)=='|');
Similarity=sum(g_alignment_P29600_P41363(2,:)==':');
g_alignment_P29600_P41363_Indentity = Indentity*100/g_alignment_P29600_P41363_length
g_alignment_P29600_P41363_Similarity = (Indentity + Similarity)*100/g_alignment_P29600_P41363_length
%% Subsection 3,4- Local alignment

disp('Subsection 3,4- Local alignment')

% Local score and alignment
[L_score_P29600_P41363,L_alignment_P29600_P41363] = swalign(sequence1,sequence2,'Alphabet','AA','ScoringMatrix','BLOSUM62')

% Alignment Length
L_alignment_P29600_P41363_length = length(L_alignment_P29600_P41363)

% Identity and Similarity Percent
Indentity=sum(L_alignment_P29600_P41363(2,:)=='|');
Similarity=sum(L_alignment_P29600_P41363(2,:)==':');
L_alignment_P29600_P41363_Indentity = Indentity*100/L_alignment_P29600_P41363_length
L_alignment_P29600_P41363_Similarity = (Indentity + Similarity)*100/L_alignment_P29600_P41363_length


%% Q5b: Alignments, scores, etc %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Q5b')

%% Sunsections 1-2, Global alignment analysis

disp('Sunsections 1-2, Global alignment analysis')

% peptidase 2-HUMAN Tripeptidyl_TPP2|P29144 
sequence3 = 'MATAATEEPFPFHGLLPKKETGAASFLCRYPEYDGRGVLIAVLDTGVDPGAPGMQVTTDGKPKIVDIIDTTGSGDVNTATEVEPKDGEIVGLSGRVLKIPASWTNPSGKYHIGIKNGYDFYPKALKERIQKERKEKIWDPVHRVALAEACRKQEEFDVANNGSSQANKLIKEELQSQVELLNSFEKKYSDPGPVYDCLVWHDGEVWRACIDSNEDGDLSKSTVLRNYKEAQEYGSFGTAEMLNYSVNIYDDGNLLSIVTSGGAHGTHVASIAAGHFPEEPERNGVAPGAQILSIKIGDTRLSTMETGTGLIRAMIEVINHKCDLVNYSYGEATHWPNSGRICEVINEAVWKHNIIYVSSAGNNGPCLSTVGCPGGTTSSVIGVGAYVSPDMMVAEYSLREKLPANQYTWSSRGPSADGALGVSISAPGGAIASVPNWTLRGTQLMNGTSMSSPNACGGIALILSGLKANNIDYTVHSVRRALENTAVKADNIEVFAQGHGIIQVDKAYDYLVQNTSFANKLGFTVTVGNNRGIYLRDPVQVAAPSDHGVGIEPVFPENTENSEKISLQLHLALTSNSSWVQCPSHLELMNQCRHINIRVDPRGLREGLHYTEVCGYDIASPNAGPLFRVPITAVIAAKVNESSHYDLAFTDVHFKPGQIRRHFIEVPEGATWAEVTVCSCSSEVSAKFVLHAVQLVKQRAYRSHEFYKFCSLPEKGTLTEAFPVLGGKAIEFCIARWWASLSDVNIDYTISFHGIVCTAPQLNIHASEGINRFDVQSSLKYEDLAPCITLKNWVQTLRPVSAKTKPLGSRDVLPNNRQLYEMVLTYNFHQPKSGEVTPSCPLLCELLYESEFDSQLWIIFDQNKRQMGSGDAYPHQYSLKLEKGDYTIRLQIRHEQISDLERLKDLPFIVSHRLSNTLSLDIHENHSFALLGKKKSSNLTLPPKYNQPFFVTSLPDDKIPKGAGPGCYLAGSLTLSKTELGKKADVIPVHYYLIPPPTKTKNGSKDKEKDSEKEKDLKEEFTEALRDLKIQWMTKLDSSDIYNELKETYPNYLPLYVARLHQLDAEKERMKRLNEIVDAANAVISHIDQTALAVYIAMKTDPRPDAATIKNDMDKQKSTLVDALCRKGCALADHLLHTQAQDGAISTDAEGKEEEGESPLDSLAETFWETTKWTDLFDNKVLTFAYKHALVNKMYGRGLKFATKLVEEKPTKENWKNCIQLMKLLGWTHCASFTENWLPIMYPPDYCVF';

% Global score and alignment
[g_score_P29600_P29144 ,g_alignment_P29600_P29144 ] = nwalign(sequence1,sequence3,'Alphabet','AA','ScoringMatrix','BLOSUM62')

% Alignment Length
g_alignment_P29600_P29144_length = length(g_alignment_P29600_P29144)

% Identity and Similarity Percent
Indentity=sum(g_alignment_P29600_P29144(2,:)=='|');
Similarity=sum(g_alignment_P29600_P29144(2,:)==':');
g_alignment_P29600_P29144_Indentity = Indentity*100/g_alignment_P29600_P29144_length
g_alignment_P29600_P29144_Similarity = (Indentity + Similarity)*100/g_alignment_P29600_P29144_length

%% Subsection 3- Local alignment

disp('Subsection 3- Local alignment')

% Local score and alignment
[L_score_P29600_P29144,L_alignment_P29600_P29144] = swalign(sequence1,sequence3,'Alphabet','AA','ScoringMatrix','BLOSUM62')

% Alignment Length
L_alignment_P29600_P29144_length = length(L_alignment_P29600_P29144)

% Identity and Similarity Percent
Indentity=sum(L_alignment_P29600_P29144(2,:)=='|');
Similarity=sum(L_alignment_P29600_P29144(2,:)==':');
L_alignment_P29600_P29144_Indentity = Indentity*100/L_alignment_P29600_P29144_length
L_alignment_P29600_P29144_Similarity = (Indentity + Similarity)*100/L_alignment_P29600_P29144_length
