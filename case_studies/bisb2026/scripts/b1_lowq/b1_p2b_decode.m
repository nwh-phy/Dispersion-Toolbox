function theta = b1_p2b_decode(u, nodes, nBranch)
%B1_P2B_DECODE  Unpack the branches' log E0^2 and log Gamma node values.
if nargin < 3, nBranch = 2; end
ne = numel(nodes.E0);
ng = numel(nodes.G);
assert(numel(u) == nBranch * (ne + ng), 'b1_p2b_decode:Size', ...
    'Expected %d branches of %d E0 and %d Gamma nodes.', nBranch, ne, ng);
theta = struct('E0sq', exp(reshape(u(1:nBranch*ne), ne, nBranch).'), ...
    'G', exp(reshape(u(nBranch*ne+1:end), ng, nBranch).'));
end
