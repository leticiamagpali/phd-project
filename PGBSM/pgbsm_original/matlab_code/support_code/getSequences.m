function SEQI = getSequences(data,TAXA,CM,Codon)

% modified for my version 16 November 2021

nL = max(setdiff(CM(:,1),CM(:,3)));

SEQI{nL} = [];
discardIDX = [];

for N = 1:nL
    
    idx1 = strfind(data,TAXA{N}) + length(TAXA{N});
    
    if N <= nL - 1
        idx2 = strfind(data,TAXA{N+1})-1;
    else
        idx2 = length(data);
    end
    
    str = data(idx1:idx2);
    str(strfind(str,newline)) = [];
    str(strfind(str,char(13))) = [];
    str(strfind(str,char(32))) = [];
    
    n_cod = length(str)/3;
    for site = 1:n_cod
        
        cID = find(strcmp(Codon,str(3*site-2:3*site)) == 1);
        
        if isempty(cID)
            SEQI{N}(site) = nan;
            discardIDX = [discardIDX;site]; %#ok<AGROW>
        else
            SEQI{N}(site) = cID;
        end
        
    end
    
    
end

% remove site patterns with missing data

 discardIDX = unique(discardIDX);
 
 for N = 1:nL
     SEQI{N}(discardIDX) = [];
 end

%% END
    
    
    
    
    
    
    
    
    
    
    
    
    
