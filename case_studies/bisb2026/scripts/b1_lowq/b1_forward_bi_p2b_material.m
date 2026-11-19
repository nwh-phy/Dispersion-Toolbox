function M = b1_forward_bi_p2b_material(S,qMin)
%B1_FORWARD_BI_P2B_MATERIAL  Source L, optionally restricted to common Q support.
% The fitted interpolation is exact on defined knots and between them.
% After restriction, the first common knot is extended as an endpoint.
c = S.cfg;
F = S.res(1).fit;
nodes = c.nodes;
nBranch = size(F.E0,1);
a = reshape(F.a,numel(nodes.A),nBranch).';
th = b1_p2b_decode(F.u,nodes,nBranch);
if nargin>1
    [th.E0sq,nodes.E0] = restrict(th.E0sq,nodes.E0,qMin,true);
    [th.G,nodes.G] = restrict(th.G,nodes.G,qMin,true);
    [a,nodes.A] = restrict(a,nodes.A,qMin,false);
end
intr = qe_lowq_joint_forward('intrinsic','Qn_E0',nodes.E0,'Qn_G',nodes.G, ...
    'Qn_A',nodes.A,'n_branch',nBranch,'ext',c.variants(1).ext);
M = struct('theta',th,'a',reshape(a.',[],1),'intrinsic',intr);
end

function [v,q] = restrict(v,q,qMin,useLog)
if qMin<=q(1), return; end
assert(qMin<q(end),'No shared momentum range.');
keep = q>qMin;
if useLog, v = log(v); end
v = [interp1(log(q),v.',log(qMin)).',v(:,keep)];
if useLog, v = exp(v); end
q = [qMin q(keep)];
end
