classdef C20V5Test < matlab.unittest.TestCase
 methods(TestClassSetup)
  function paths(t)
   root=fileparts(fileparts(mfilename('fullpath')));
   t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src'),IncludingSubfolders=true));
   t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'case_studies','bisb2026','scripts')));
  end
 end
 methods(Test)
  function areaScopes(t)
   E=(300:4:2100)'; Y=exp(-E/1000);
   a=qe_area_scopes(E(E<=1800),Y(E<=1800)); b=qe_area_scopes(E(E<=2000),Y(E<=2000)); c=qe_area_scopes(E,Y);
   t.verifyEqual(a.area_reference_window,b.area_reference_window); t.verifyEqual(a.area_reference_window,c.area_reference_window);
   t.verifyGreaterThan(c.area_fit_window,b.area_fit_window); t.verifyGreaterThan(b.area_fit_window,a.area_fit_window);
  end
  function noExtrapolation(t)
   E=(400:4:1600)'; a=qe_area_scopes(E,ones(size(E)));
   t.verifyFalse(a.reference_window_fully_observed); t.verifyEqual(a.area_reference_window,1200);
   t.verifyEqual(a.reference_actual_support_meV,[400 1600]);
  end
  function zeroSlopeAndJacobian(t)
   E=(300:8:1800)'; q=.001:.0005:.003; p=[1.2 .1*ones(1,5) .9 .9 .4 .7*ones(1,5)];
   [Y,~,~,J]=qe_member_prediction(E,q,p,1,'lorentz_symmetric',1000);
   pp=p; pp(7)=pp(7)+1e-6; y2=qe_member_prediction(E,q,pp,1,'lorentz_symmetric',1000);
   pp(7)=p(7)-1e-6; ym=qe_member_prediction(E,q,pp,1,'lorentz_symmetric',1000);
   t.verifyEqual(Y(:,1),Y(:,end)); t.verifyEqual((y2(:)-ym(:))/2e-6,J(:,7),AbsTol=2e-5);
  end
  function dlJacobian(t)
   E=(300:8:1800)'; q=.001:.0005:.003; p=[1.2 .1*ones(1,5) .8 1.1 .4 .7*ones(1,5)];
   [Y,~,~,J]=qe_member_prediction(E,q,p,1,'lorentz',1e6);
   pp=p; pp(9)=pp(9)+1e-6; y2=qe_member_prediction(E,q,pp,1,'lorentz',1e6);
   t.verifyEqual((y2(:)-Y(:))/1e-6,J(:,9),AbsTol=3e-5);
  end
  function recoverLocalSlope(t)
   E=(300:12:1800)'; q=.001:.0005:.003; p=[1.2 .1*ones(1,5) .8 1.0 .4 .7*ones(1,5)];
   Y=qe_member_prediction(E,q,p,1,'lorentz_symmetric',1000);
   f=qe_fit_member_models(E,q,Y,n_starts=4);
   t.verifyTrue(all([f.success])); t.verifyLessThan(f(1).objective,1e-7);
   t.verifyEqual(f(1).slope_meV_A,100000,AbsTol=100);
   t.verifyEqual(f(2).witness.prediction,f(1).prediction,AbsTol=1e-9);
   t.verifyEqual(f(1).observation_count,numel(Y)); t.verifyEqual(numel(f(2).candidates),4);
  end
  function scaleAndEffectiveOptions(t)
   E=(300:12:1800)'; q=.001:.0005:.003; p=[1.2 .1*ones(1,5) .8 1.0 .4 .7*ones(1,5)];
   Y=qe_member_prediction(E,q,p,1,'lorentz_symmetric',1000);
   a=qe_fit_member_models(E,q,Y,n_starts=2,energy_window=[400 1600]); b=qe_fit_member_models(E,q,Y*3,n_starts=2,energy_window=[400 1600]);
   t.verifyEqual(a(1).effective_options.energy_window,[400 1600]); t.verifyEqual(numel(a(1).candidates),2);
   t.verifyEqual(b(1).prediction,a(1).prediction*3,AbsTol=1e-7);
  end
 end
end
