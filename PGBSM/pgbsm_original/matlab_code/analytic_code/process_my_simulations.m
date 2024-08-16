% batch-process alignments simulated
% by C T Jones August 2022

clc
clear
close all

global path_to_pgbsm

path_to_pgbsm = pwd;% sets path_to_pgbsm

path_to_data = [path_to_pgbsm '\simulated_alignments\'];

% make paths so Matlab can find things

addpath(path_to_data)
addpath([path_to_pgbsm '\genetic_codes'])
addpath([path_to_pgbsm '\matlab_code\support_code'])
addpath([path_to_pgbsm '\matlab_code\likelihood_functions'])

% load Codon and AminoAcid data structures

genetic_code = 'Mammalian_Mitochondrial_GeneticCode';
load(genetic_code);

% get common data

fileList_tree = dir([path_to_data '\tree_data*']);
fileList_phenotypeMap = dir([path_to_data '\phenotype_map*']);
fileList_labels = dir([path_to_data '\taxon_labels*']);

[T,nL,CM] = readTreeFile(path_to_data,fileList_tree.name);
[Y,pi_state] = readPhenotypeMap(path_to_data,fileList_phenotypeMap.name);
Labels = getLabels(path_to_data,fileList_labels.name);

% get a list of files to process

fileList = dir([path_to_data '\seqfile*.txt']);
Nfiles = size(fileList,1);

%% fit models

% RaMoSS
Sim_Output(Nfiles).nulRaMoSS.mle = mle0;
Sim_Output(Nfiles).nulRaMoSS.LL = -LL0;
Sim_Output(Nfiles).altRaMoSS.mle = mle1;
Sim_Output(Nfiles).altRaMoSS.LL = -LL1;
Sim_Output(Nfiles).altRaMoSS.Post = PostH;

% Null M0
Sim_Output(Nfiles).Nul.mle = null_mle;
Sim_Output(Nfiles).Nul.LL = -LL;

% PG-BSM BW
Sim_Output(Nfiles).BW.mle = mle;
Sim_Output(Nfiles).BW.LL = -LL;
Sim_Output(Nfiles).BW.POST = POST;

% PG-BSM CW
Sim_Output(Nfiles).CW.mle = mle;
Sim_Output(Nfiles).CW.LL = -LL;
Sim_Output(Nfiles).CW.POST = POST;

% PG-BSM rCW
Sim_Output(Nfiles).rCW.mle = mle;
Sim_Output(Nfiles).rCW.LL = -LL;
Sim_Output(Nfiles).rCW.POST = POST;

for file = 1:Nfiles
    
    % construct alignment data structure
    data = readSequenceFile(path_to_data,fileList(file).name);
    SEQI = getSequences(data,Labels,CM,Codon);
    
    % construct matrix of target frequencies using F3x4
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
    Sim_Output(file).nulRaMoSS.mle = mle0;
    Sim_Output(file).nulRaMoSS.LL = -LL0;
    Sim_Output(file).altRaMoSS.mle = mle1;
    Sim_Output(file).altRaMoSS.LL = -LL1;
    
    save([path_to_pgbsm '\simulated_alignments\Sim_Output'],'Sim_Output')
    
    % ***********************************************************
    % compute posterior probability of being a heterotachous site
    % added for my own version of the code, 15 September 2021
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
    
    Sim_Output(file).altRaMoSS.Post = PostH;
    
    save([path_to_pgbsm '\simulated_alignments\Sim_Output'],'Sim_Output')
    
    %% fit null PG-BSM, mle = (piM3 w1CL w2CL p1CL delta kappa lambda B)
    
    [null_mle,LL] = fitNullPGBSM(CM,SEQI,TF,Y,pi_state,genetic_code);
    
    Sim_Output(file).Nul.mle = null_mle;
    Sim_Output(file).Nul.LL = -LL;
    
    save([path_to_pgbsm '\simulated_alignments\Sim_Output'],'Sim_Output')
    
    lambda = null_mle(5);
    CM(:,2) = null_mle(8:end);
    
    %% fit alternative PG-BSM BW
    
    [mle,LL,POST] = fitPGBSM_BW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
    Sim_Output(file).BW.mle = mle;
    Sim_Output(file).BW.LL = -LL;
    Sim_Output(file).BW.POST = POST;
    
    save([path_to_pgbsm '\simulated_alignments\Sim_Output'],'Sim_Output')
    
    %% fit alternative PG-BSM CW
    
    [mle,LL,POST] = fitPGBSM_CW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
    Sim_Output(file).CW.mle = mle;
    Sim_Output(file).CW.LL = -LL;
    Sim_Output(file).CW.POST = POST;
    
    save([path_to_pgbsm '\simulated_alignments\Sim_Output'],'Sim_Output')
    
    %% fit alternative PG-BSM rCW
    
    [mle,LL,POST] = fitPGBSM_rCW(null_mle,CM,SEQI,TF,Y,pi_state,lambda,genetic_code);
    Sim_Output(file).rCW.mle = mle;
    Sim_Output(file).rCW.LL = -LL;
    Sim_Output(file).rCW.POST = POST;
    
    save([path_to_pgbsm '\simulated_alignments\Sim_Output'],'Sim_Output')
    
end

disp('Done!')

%% END