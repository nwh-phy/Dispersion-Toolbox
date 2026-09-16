classdef C20PhysicsPreviewTest < matlab.unittest.TestCase
 methods(TestClassSetup)
  function paths(t)
   t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'),IncludingSubfolders=true));
  end
 end
 methods(Test)
  function centroidIdentity(t)
   E=(300:4:1800)'; m=peak_models('lorentz_symmetric'); C=[m.model_fn(600,800,1000,E),m.model_fn(1200,1000,800,E)];
   r=qe_spectral_centroids(E,C,C); t.verifyLessThan(r.identity_error,1e-10); t.verifyEqual(r.frozen_centroid,r.combined_centroid,AbsTol=1e-9);
  end
  function reweightWithoutShapeChange(t)
   E=(300:4:1800)'; m=peak_models('lorentz_symmetric'); C=[m.model_fn(600,800,1000,E),m.model_fn(1200,1000,800,E)];
   a=qe_spectral_centroids(E,C,C); b=qe_spectral_centroids(E,C.*[.1 1],C);
   t.verifyGreaterThan(b.frozen_centroid,a.frozen_centroid); t.verifyEqual(b.combined_centroid,b.frozen_centroid,AbsTol=1e-9);
  end
 end
end
