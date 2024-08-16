% process user data
% by C T Jones August 2022

clc
clear
close all

global path_to_pgbsm

path_to_pgbsm = pwd;% sets path_to_pgbsm

% make paths so Matlab can find things
addpath([path_to_pgbsm '\genetic_codes'])
addpath([path_to_pgbsm '\matlab_code\support_code'])
addpath([path_to_pgbsm '\matlab_code\likelihood_functions'])

% path to data
path_to_data = [path_to_pgbsm '\usersReal_alignments\'];
addpath(path_to_data)

%% BEGIN USER INPUT

% choose a genetic code (cut and paste)
%
% Standard_GeneticCode
% Mammalian_Mitochondrial_GeneticCode
% Invertebrate_Mitochondrial_GeneticCode

genetic_code = 'Mammalian_Mitochondrial_GeneticCode';
load(genetic_code);

% specify file names
tree = 'formatted_tree.txt';
sequence = 'formatted_sequences.txt';
phenotype = 'phenotype_map.txt';
labels = 'taxon_labels.txt';

%% END USER INPUT

%% read in data

[T,nL,CM] = readTreeFile(path_to_data,tree);             % tree structure
[Y,pi_state] = readPhenotypeMap(path_to_data,phenotype); % phenotype map

TAXA = getLabels(path_to_data,labels);                   % taxa labels

data = readSequenceFile(path_to_data,sequence);          % sequence data
SEQI = getSequences(data,TAXA,CM,Codon);

%% construct matrix of target frequencies using F3x4

pi_nuc = getNucFreq(SEQI,Codon);
TF = makeTargetFrequencies(pi_nuc,Codon);

%% fit to M0 to obtain an initial estimate of branch lengths

CM(:,2) = ones(size(CM(:,2)));

% mle = (w, kappa, B)
[mle,LLM0] = fitM0(CM,SEQI,TF,AminoAcid,Codon);

Branch_Lengths = mle(3:end);
CM(:,2) = Branch_Lengths(:);

%% fit null and alternate RaMoSS

% mle = (propCL w1M3 w2M3 p1M3 w1CL w2CL p1CL delta kappa B)

[mle0,mle1,LL0,LL1] = fitRaMoSS(Branch_Lengths,CM,SEQI,TF,genetic_code);
Output.nulRaMoSS.mle = mle0;
Output.nulRaMoSS.LL = -LL0;
Output.altRaMoSS.mle = mle1;
Output.altRaMoSS.LL = -LL1;

save([path_to_pgbsm '\usersReal_alignments\Output'],'Output')

% ***********************************************************
% compute posterior probability of being a heterotachous site
% added by C T Jones on 15 September 2021
% **********************************************************

piCL = mle1(1);
n_sites = length(SEQI{1});

PostH = zeros(1,n_sites);

disp('Computing RaMoSS posteriors ...')
for site = 1:n_sites
    
    tempSEQI = cell(1,nL);
    
    for n = 1:nL
        tempSEQI{n} = SEQI{n}(site);
    end
    
    [~,~,LM3,LCL] = likelihoodRaMoSS(mle1,CM,TF,tempSEQI,genetic_code,0);
    PostH(site) = LCL*piCL/(LM3*(1-piCL) + LCL*piCL);
    
end

Output.altRaMoSS.Post = PostH;

save([path_to_pgbsm '\usersReal_alignments\Output'],'Output')

%% fit null PG-BSM, mle = (piM3 w1CL w2CL p1CL delta kappa lambda B)

[null_mle,LL] = fitNullPGBSM(CM,SEQI,TF,Y,pi_state,genetic_code);

Output.Nul.mle = null_mle;
Output.Nul.LL = -LL;

save([path_to_pgbsm '\usersReal_alignments\Output'],'Output')

lambda = null_mle(5);
CM(:,2) = null_mle(8:end);

%% fit alternative PG-BSM BW

[mle,LL,POST] = fitPGBSM_BW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
Output.BW.mle = mle;
Output.BW.LL = -LL;
Output.BW.POST = POST;

save([path_to_pgbsm '\usersReal_alignments\Output'],'Output')

%% fit alternative PG-BSM CW

[mle,LL,POST] = fitPGBSM_CW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
Output.CW.mle = mle;
Output.CW.LL = -LL;
Output.CW.POST = POST;

save([path_to_pgbsm '\usersReal_alignments\Output'],'Output')

%% fit alternative PG-BSM rCW

[mle,LL,POST] = fitPGBSM_rCW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
Output.rCW.mle = mle;
Output.rCW.LL = -LL;
Output.rCW.POST = POST;

save([path_to_pgbsm '\usersReal_alignments\Output'],'Output')

disp('Done!')

%% END