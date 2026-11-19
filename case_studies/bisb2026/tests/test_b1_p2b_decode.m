function tests = test_b1_p2b_decode
tests = functiontests(localfunctions);
end

function testArbitraryNodeCountsAndBranchOrder(t)
nodes = struct('E0',[1 2 3 4 5 6]*1e-3,'G',[1 3 5 6]*1e-3);
E0sq = [1:6; 7:12].^2;
G = [20:23; 30:33];
u = [reshape(log(E0sq).',1,[]),reshape(log(G).',1,[])];
th = b1_p2b_decode(u,nodes);
verifyEqual(t,th.E0sq,E0sq,'RelTol',1e-14);
verifyEqual(t,th.G,G,'RelTol',1e-14);
verifyError(t,@() b1_p2b_decode(u(1:end-1),nodes),'b1_p2b_decode:Size');
end
