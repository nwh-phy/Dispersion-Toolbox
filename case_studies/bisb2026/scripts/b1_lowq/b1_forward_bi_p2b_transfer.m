function X = b1_forward_bi_p2b_transfer(P,F,c)
%B1_FORWARD_BI_P2B_TRANSFER  Scale donor material on recipient training data.
% P is a fixed-material replay through the recipient operator. F is the
% recipient's own fit. Its nuisance tail, including the held block, is
% held at the training-derived values used for its own prediction.
[~,ni,nb] = size(P.Y);
tail = zeros(size(P.E));
m = P.E>=c.holdout(1);
if F.tailOn, tail(m) = (P.E(m)/c.tail.E_ref).^(-c.tail.p); end
nuisance = tail.*reshape(F.b,1,ni,nb);
train = logical(F.mask);
held = P.mask & ~train;
assert(isequal(train|held,P.mask) && ~any(train(:,:,c.lobo_block),'all') && ...
    ~any(held(:,:,setdiff(1:nb,c.lobo_block)),'all'), ...
    'b1_forward_bi_p2b_transfer:Split','Recipient masks must use the same block split.');
u = P.mu(train)./sqrt(P.S2(train));
y = (P.Y(train)-nuisance(train))./sqrt(P.S2(train));
assert(all(isfinite(u)) && dot(u,u)>0,'b1_forward_bi_p2b_transfer:Empty','No donor signal in training.');
scale = max(0,dot(u,y)/dot(u,u));
mu = scale*P.mu+nuisance;
z = (P.Y-mu)./sqrt(P.S2);
own = (P.Y-F.mu)./sqrt(P.S2);
trainN = nnz(train);
heldN = nnz(held);
crossChi2 = [sum(z(train).^2);sum(z(held).^2)];
ownChi2 = [sum(own(train).^2);sum(own(held).^2)];
counts = [trainN;heldN];
ownDof = [F.dof;heldN];
crossDof = [trainN-1;heldN];
X = struct('scale',scale,'table',table(["training";"held block"],counts, ...
    ownChi2,crossChi2,ownDof,crossDof,ownChi2./ownDof,crossChi2./crossDof, ...
    (crossChi2-ownChi2)./counts,'VariableNames',{'split','n','own_chi2','cross_chi2', ...
    'own_dof','cross_dof','own_chi2dof','cross_chi2dof','delta_chi2_per_point'}), ...
    'mu',mu,'z',z,'own_z',own);
end
