% The PG-BSM code requires a tree with taxa indicated by numbers 1 to nL 
% where nL is the number of taxa. Sequences are assumed to occur in the same 
% order in the sequence file as they appear in the tree. This script was 
% written to format a sequence file to meet this requirement.

% Taxon labels in the tree must match those labels in the sequence.

% by C T Jones August 2022

clc
clear
close all

global path_to_pgbsm

path_to_pgbsm = pwd;% sets path_to_pgbsm
addpath([path_to_pgbsm '\matlab_code\support_code'])

% make paths so Matlab can find things
addpath([path_to_pgbsm '\genetic_codes'])

% specify the path to your data
path_to_data = [path_to_pgbsm '\usersReal_alignments\'];
addpath(path_to_data)

fid = fopen([path_to_data '\my_sequences.txt']); % <- name your sequence file
data = char(fread(fid)');
fclose(fid); 

fid = fopen([path_to_data '\my_tree.txt']);% <- name your tree file
T = char(fread(fid)');
T(strfind(T,newline)) = [];
T(strfind(T,char(13))) = [];
fclose(fid); 

%% read nL and ncod

% the first line is assumed to give the number of taxa and sequence length

idx = strfind(data,newline);
str = data(1:idx(1));

data = data(idx(1):end);

idx = strfind(str,char(32));
nL = str2double(str(1:idx-1));   % taxa
nN = str2double(str(idx+1:end)); % nucleotides

%% read taxa from the tree

stop_idx = strfind(T,':');

TAXA = cell(nL,1);
counter = 1;
for n = 1:length(stop_idx)
    if ~strcmp(T(stop_idx(n)-1),')')
        str = T(1:stop_idx(n)-1);
        start_idx = union(strfind(T(1:stop_idx(n)-1),'('),strfind(T(1:stop_idx(n)-1),','));
        TAXA{counter} = T(start_idx(end)+1:stop_idx(n)-1);
        counter = counter + 1;
    end
end

fid = fopen([path_to_data 'taxon_labels.txt'],'wt+');
for n= 1:nL
    fprintf(fid,[TAXA{n} newline]);
end
fclose('all');

%% replace taxa labels with numbers

newT = T;
for n = 1:nL
    idx = strfind(newT,TAXA{n});
    newT = [newT(1:idx-1) num2str(n) newT(idx + length(TAXA{n}):end)];
end

fid = fopen([path_to_data 'formatted_tree.txt'],'wt+');
fprintf(fid,newT);
fclose('all');

%% shuffle sequences

% each sequence is assumed to start with a taxa name and end with a newline

EOL = strfind(data,newline);

fid = fopen([path_to_data 'formatted_sequences.txt'],'wt+');
fprintf(fid,[num2str(nL) char(32) num2str(nN) newline newline]);

for n = 1:nL
    
    disp(TAXA{n});
    
    idx1 = strfind(data,TAXA{n});
    idx2 = EOL(EOL > idx1(1));
    str = data(idx1:idx2(1));
    
    fprintf(fid,[str newline]);
 
end

fprintf(fid,newline);
fclose('all');

%% END