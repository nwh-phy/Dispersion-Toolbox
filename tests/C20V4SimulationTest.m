classdef C20V4SimulationTest < matlab.unittest.TestCase
 methods(TestClassSetup)
  function paths(t)
   root=fileparts(fileparts(mfilename('fullpath')));
   t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src'),IncludingSubfolders=true));
   t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'case_studies','bisb2026','scripts')));
  end
 end
 methods(Test)
  function h0DoesNotUseExperimentalMean(t)
   E=(300:4:1800)'; q=[.002 .0025 .003];
   a=c20_v4_simulate(E,q,'static_single',1); b=c20_v4_simulate(E,q,'static_single',2);
   t.verifyEqual(a.mean,b.mean); t.verifyNotEqual(a.noise,b.noise);
   t.verifyEqual(size(a.truth_parameters,1),1);
   t.verifyEqual(a.observed-a.mean,a.noise,AbsTol=1e-12);
  end
  function realMembersAndDoubleTruth(t)
   E=(300:4:1800)'; q=[.002 .0025 .003];
   a=c20_v4_simulate(E,q,'q_mixed_single',1); b=c20_v4_simulate(E,q,'overlapping_double',2);
   t.verifyEqual(a.mean,mean(a.member_truth,2)); t.verifyEqual(a.q_members,q);
   t.verifyEqual(size(b.truth_parameters,1),2);
   t.verifyTrue(contains(a.role,'not_experimental'));
  end
 end
end
