function c20_v5_controls(out)
ref=jsondecode(fileread(fullfile(out,'parent_packet_reference.json'))); parent=ref.parent;
a=load(fullfile(parent,'validation','simulation_phase1_background_mismatch.mat')); background_legacy=a.trials;
member=load(fullfile(out,'member_models','R1_lorentz_symmetric','member_fit.mat')); realM1=member.d.fits(1); assert(realM1.success);
cfg=struct('trials',10,'seed_base',609120,'member_starts',4,'background_starts',4,'background_extra',4, ...
 'noise','hypothetical iid Gaussian on members; no experimental false-positive rate','data_constrained_truth','R1 A0 symmetric M1 selected member predictions');
code_paths={mfilename('fullpath'),which('qe_fit_member_models'),which('qe_member_prediction'),which('qe_compare_component_models'),which('peak_models')}; code_paths{1}=[code_paths{1} '.m'];
hashes=cellfun(@(p)c20_v4_io('hash',p),code_paths,'UniformOutput',false); generator_hash=c20_v4_io('digest',strjoin(hashes,'|'));
cfg.generator_hash=generator_hash; cfg.data_hash=c20_v4_io('hash',fullfile(out,'member_models','R1_lorentz_symmetric','member_fit.mat'));
c20_v4_io('json',fullfile(out,'validation','control_config.json'),cfg); rows={};
for scenario=1:3
 trials=cell(1,cfg.trials);
 for k=1:cfg.trials
  seed=cfg.seed_base+100*scenario+k; stream=RandStream('mt19937ar','Seed',seed);
  if scenario==1
   old=background_legacy{k}; s=old.simulation;
   fits=qe_compare_component_models(s.E,s.observed,n_starts=cfg.background_starts,extra_starts=cfg.background_extra, ...
    start_policy='independent',seed=40912,energy_window=[300 1800],baseline_mode='power_law_plus_nonnegative_constant');
   improvement=(fits(1).normalized_sse-fits(2).normalized_sse)/fits(1).normalized_sse;
   areas=trapz(s.E,fits(2).components); flag=improvement>.1&&min(areas)/sum(areas)>.05&&diff(fits(2).parameters(:,1))>4;
   err=sqrt(mean((fits(1).prediction-s.mean).^2));
   trials{k}=struct('observed',s.observed,'E',s.E,'truth',s.mean,'fits',fits,'legacy_B0_fits',old.fits,'seed',s.seed, ...
    'generator_hash',generator_hash,'data_hash',cfg.data_hash);
   rows(end+1,:)={string('background_correct_B1'),k,all([fits.success]),improvement,flag,err}; %#ok<AGROW>
  else
   if scenario==2
    E=realM1.E; q=[.002 .0025 .003]; pm=[1.3 .12*ones(1,3) .74 1.06 .28 .9*ones(1,3)];
    truth=qe_member_prediction(E,q,pm,1,'lorentz_symmetric',1000); name='legacy_qmix_member_positive';
   else
    E=realM1.E; q=realM1.q; truth=realM1.prediction; name='R1_data_constrained_M1';
   end
   sigma=.01*max(truth,[],'all'); noise=sigma*randn(stream,size(truth)); observed=truth+noise;
   fits=qe_fit_member_models(E,q,observed,n_starts=cfg.member_starts,seed=40912,energy_window=[300 1800]);
   improvement=(fits(1).objective-fits(2).objective)/fits(1).objective;
   err=sqrt(mean((fits(1).prediction-truth).^2,'all'))/max(truth,[],'all');
   trials{k}=struct('observed',observed,'truth',truth,'noise',noise,'sigma',sigma,'q',q,'E',E,'fits',fits,'seed',seed, ...
    'generator_hash',generator_hash,'data_hash',cfg.data_hash);
   rows(end+1,:)={string(name),k,all([fits.success]),improvement,NaN,err}; %#ok<AGROW>
  end
  disp(sprintf('CONTROL scenario%d %d/%d',scenario,k,cfg.trials));
 end
 save(fullfile(out,'validation',sprintf('control_scenario%d.mat',scenario)),'trials','cfg','-v7');
end
writetable(cell2table(rows,'VariableNames',{'scenario','trial','success','relative_objective_gain','legacy_engineering_split_not_detection','M1_truth_curve_RMS'}),fullfile(out,'validation','positive_controls.csv'));
end
