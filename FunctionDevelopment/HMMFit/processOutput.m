basepaths_path = '../basepaths.mat';
load(basepaths_path)
for ii = 1:length(basepaths)
    basepath = basepaths{ii};
    fprintf('%d/%d\n', ii, length(basepaths))
    if exist( fullfile( basepath, 'GammaProcessed' ) )
        
        cd(fullfile( basepath, 'GammaProcessed' ))
        fils = dir(pwd); fils(1:2) = [];
        
        NREMall = [];
        WAKEall = [];
        for kp = 1:length(fils)
            
            if fils(kp).bytes == 128
                delete(fils(kp).name)
            else
                v = load(fils(kp).name);
                % Store NREM
                if isfield(v, 'NREM')
                    NREMall = [NREMall v.NREM];
                end
                % Store WAKE
                if isfield(v, 'WAKE')
                    WAKEall = [WAKEall v.WAKE];
                end
            end
            
        end
        save('hmm_out.mat', 'WAKEall', 'NREMall')
    end
    
end