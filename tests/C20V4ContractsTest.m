classdef C20V4ContractsTest < matlab.unittest.TestCase
 methods(TestClassSetup)
  function paths(t)
   root=fileparts(fileparts(mfilename('fullpath')));
   t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src'),IncludingSubfolders=true));
   t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'case_studies','bisb2026','scripts')));
  end
 end
 methods(Test)
  function rawHighMapsToP2(t)
   p=[0 2 1.3 5 1 .7 .2 2];
   m=qe_component_mapping(p,[0 0 .3 .004 0 .3 .004 0],[Inf 6 1.8 5 Inf 1.8 5 Inf],2,10,1000,'lorentz_symmetric');
   t.verifyEqual(m.order,[2;1]); t.verifyEqual(m.ordered_component(4),2);
   t.verifyTrue(m.upper_hit(4)); t.verifyEqual(m.parameters(2,2),5000);
   t.verifyEqual(m.ordered_component(1),0); t.verifyTrue(m.lower_hit(1));
  end
  function nanAndEqualCenter(t)
   p=[0 0 .7 .2 0 .7 .3 1];
   m=qe_component_mapping(p,zeros(1,8),inf(1,8),2,1,1000,'lorentz_symmetric');
   bad=qe_component_mapping(nan(1,5),zeros(1,5),inf(1,5),1,1,1000,'lorentz_symmetric');
   t.verifyTrue(m.equal_center); t.verifyTrue(m.lower_hit(5)); t.verifyFalse(bad.finite);
   t.verifyTrue(all(bad.boundary_type=="not_assessable"));
  end
  function centeredGeometry(t)
   q=C20V4ContractsTest.fixture(); b=qe_centered_bins(q,.0025,[1 3 5]);
   t.verifyEqual([b.q_Ainv],repmat(.0025,1,3),AbsTol=1e-12);
   t.verifyEqual(b(3).source_channel,103:107); t.verifyEqual(b(3).q_width,.0025,AbsTol=1e-12);
   t.verifyEqual(b(3).sum,5*b(3).mean); t.verifyTrue(all(isnan(b(3).variance_mean)));
   t.verifyFalse(b(3).independent_of_other_N);
  end
  function rejectsGapAndMask(t)
   q=C20V4ContractsTest.fixture(); q.intensity(2,4)=NaN;
   b=qe_centered_bins(q,.0025,3); t.verifyFalse(b.valid);
   q=C20V4ContractsTest.fixture(); q.source_channel(5:end)=q.source_channel(5:end)+1;
   b=qe_centered_bins(q,.0025,3); t.verifyFalse(b.valid);
  end
  function rejectsCentralAndIrregular(t)
   q=C20V4ContractsTest.fixture(); q.q_Ainv=q.q_Ainv-.0025;
   b=qe_centered_bins(q,0,3); t.verifyFalse(b.valid);
   q=C20V4ContractsTest.fixture(); q.q_Ainv(4)=q.q_Ainv(4)+.0001;
   b=qe_centered_bins(q,.0025,3); t.verifyFalse(b.valid);
  end
  function noncircularSignsAndSupport(t)
   E=(-20:20)'; X=zeros(41,3); X(18,1)=10; X(21,2)=10; X(25,3)=10; X(1,:)=1; X(end,:)=2;
   a=qe_zlp_integer_align(E,X,X,[-10 10]);
   [~,k]=max(a.aligned,[],1);
   t.verifyEqual(a.measured_offset_pixels,[-3 0 4]); t.verifyEqual(a.correction_pixels,[3 0 -4]);
   t.verifyEqual(k,repmat(find(a.E==0),1,3));
   t.verifyEqual(sum(a.aligned,1),[sum(X(1:34,1)),sum(X(4:37,2)),sum(X(8:41,3))]);
   t.verifyEqual(a.support,(4:37)');
  end
  function invalidZlpIsUnknown(t)
   E=(-20:20)'; X=zeros(41,2); X(21,1)=10;
   a=qe_zlp_integer_align(E,X,X,[-10 10]);
   t.verifyFalse(a.valid(2)); t.verifyTrue(isnan(a.correction_pixels(2)));
   t.verifyTrue(all(isnan(a.aligned(:,2))));
   t.verifyError(@()qe_zlp_integer_align(E,zeros(41,2),X,[-10 10]),'qe_zlp_integer_align:NoZLP');
  end
  function nestedAndIndependentStreams(t)
   E=(300:16:1800)'; model=peak_models('lorentz_symmetric'); Y=.1*(E/1000).^-1.3+model.model_fn(850,300,400,E);
   a=qe_compare_component_models(E,Y,n_starts=3,start_policy='independent');
   b=qe_compare_component_models(E,Y,n_starts=3,h0_n_starts=5,start_policy='independent');
   t.verifyEqual(vertcat(a(2).candidates.p0),vertcat(b(2).candidates.p0));
   t.verifyTrue(a(2).witness.available); t.verifyEqual(a(2).witness.prediction,a(1).prediction,AbsTol=1e-10);
   t.verifyEqual(a(2).witness.candidate_type,'feasible_zero_amplitude_not_optimized');
  end
  function backgroundAndScale(t)
   E=(300:16:1800)'; p=[.1 1.2 .8 .2 .4];
   y=qe_component_prediction(E,p,'lorentz_symmetric',1,1000);
   yc=qe_component_prediction(E,[p .2],'lorentz_symmetric',1,1000);
   t.verifyEqual(yc-y,.2*ones(size(E)),AbsTol=1e-14);
   a=qe_compare_component_models(E,yc,n_starts=3,baseline_mode='power_law_plus_nonnegative_constant');
   b=qe_compare_component_models(E,3*yc,n_starts=3,baseline_mode='power_law_plus_nonnegative_constant');
   t.verifyEqual(numel(a(1).lb),6); t.verifyEqual(numel(a(2).lb),9);
   t.verifyEqual(b(1).prediction,3*a(1).prediction,AbsTol=1e-5);
   t.verifyEqual(a(2).witness.prediction,a(1).prediction,AbsTol=1e-10);
  end
  function failedCandidatesExport(t)
   E=(300:16:1800)'; z=qe_compare_component_models(E,exp(-E/500),n_starts=1,max_iterations=0);
   t.verifyFalse(any([z.success])); t.verifyTrue(all(isnan([z.selected_start])));
   folder=t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
   u=struct('bin_size_requested',3,'source_q_count',3,'q_Ainv',.0025,'q_left',.00175,'q_right',.00325);
   stats=c20_v4_export_fits({struct('key','failure_test','unit',u,'fits',z)},folder.Folder);
   t.verifyEqual(stats.failed_models,2); t.verifyEqual(stats.components,3);
  end
 end
 methods(Static)
  function q=fixture()
   q=struct('q_Ainv',.0005:.0005:.0045,'dq_Ainv',.0005,'source_channel',101:109, ...
    'energy_meV',(300:4:1800)','intensity',ones(376,9));
  end
 end
end
