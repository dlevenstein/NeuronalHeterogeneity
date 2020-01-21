function dist=KLDiv(P,Q,varargin)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
%
%Options
%   'symmetric'     true/false (default: false)
%   'epsilon'       small amount added to 0 bins (default 1e-8)
%%
p = inputParser;
addParameter(p,'epsilon',1e-8)
addParameter(p,'symmetric',false)


parse(p,varargin{:})
epsilon = p.Results.epsilon;
symmetric = p.Results.symmetric;

%%
P(P<=0) = epsilon;
Q(Q<=0) = epsilon;

if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

if sum(isinf(P(:))) + sum(isinf(Q(:)))
   error('the inputs contain non-finite values!') 
end


%Pass through NaNs
passnans = any(isnan(P),2);
if size(Q,1)==1
    if any(isnan(Q))
        dist = nan(size(passnans));
        return
    end
elseif size(Q,1)==size(P,1)
    passnans = passnans | any(isnan(Q),2);
end

% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
    if symmetric
        temp =  repmat(Q,[size(P,1) 1]).*log(repmat(Q,[size(P,1) 1])./P);
        temp(isnan(temp))=0;% resolving the case when P(i)==0
        dist = dist+sum(temp,2);
    end
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = sum(temp,2);
    
    if symmetric
    temp =  Q.*log(Q./P);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = dist+sum(temp,2);
    end
    
end

dist(passnans) = nan;
