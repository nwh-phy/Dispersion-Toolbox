function theta = b1_p2b_decode(u, nodes)
%B1_P2B_DECODE  Unpack the two branches' log E0^2 and log Gamma node values.
ne = numel(nodes.E0);
ng = numel(nodes.G);
assert(numel(u) == 2 * (ne + ng), 'b1_p2b_decode:Size', ...
    'Expected two branches of %d E0 and %d Gamma nodes.', ne, ng);
theta = struct('E0sq', exp(reshape(u(1:2*ne), ne, 2).'), ...
    'G', exp(reshape(u(2*ne+1:end), ng, 2).'));
end
