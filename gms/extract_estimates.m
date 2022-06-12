function [X,N,L,N_c]=extract_estimates(glmb,model)
%extract estimates via best cardinality, then 
%best component/hypothesis given best cardinality, then
%best means of tracks given best component/hypothesis and cardinality
[~,mode] = max(glmb.cdn);
N = mode-1;
X= zeros(model.x_dim,N);
L= zeros(2,N);
[~,idxcmp]= max(glmb.w.*(glmb.n==N));
for n=1:N
    [~,idxtrk]= max(glmb.tt{glmb.I{idxcmp}(n)}.w);
    X(:,n)= glmb.tt{glmb.I{idxcmp}(n)}.m(:,idxtrk);
    L(:,n)= glmb.tt{glmb.I{idxcmp}(n)}.l;
end
N_c = glmb.clutter(idxcmp) ; 
end