function parms = EnterDefaultParms(parms,defaultParms)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%DLevenstein 2019
%%
if isempty(parms)
    parms = defaultParms;
end
parmfields = fieldnames(defaultParms);

for pp = 1:length(parmfields)
    %If the field isn't in the paramters strucutre, fill it with the default
    if ~isfield(parms,parmfields{pp}) 
        parms.(parmfields{pp}) = defaultParms.(parmfields{pp});
    end
end

end

