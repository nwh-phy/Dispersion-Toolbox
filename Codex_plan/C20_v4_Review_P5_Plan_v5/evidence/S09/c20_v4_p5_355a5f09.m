function c20_v4_p5(out)
% Four fixed scenarios: 10 smoke + 100 independent pilot draws each.
dest=fullfile(out,'validation'); assert(isfolder(dest));
src={mfilename('fullpath'),which('c20_v4_simulate'),which('c20_v4_profiles'),which('C20V4SimulationTest')};
src{1}=[src{1} '.m']; src{4}=fullfile(pwd,'tests','C20V4SimulationTest.m');
provenance=cell(numel(src),2);
for k=1:numel(src)
 provenance(k,:)={string(src{k}),string(c20_v4_io('hash',src{k}))};
 [~,name,ext]=fileparts(src{k}); h=char(provenance{k,2}); copyfile(src{k},fullfile(out,'provenance','source_snapshot',[name '_' h(1:8) ext]));
end
writetable(cell2table(provenance,'VariableNames',{'path','sha256'}),fullfile(out,'provenance',['P5_stage_source_hashes_' char(datetime('now','Format','HHmmssSSS')) '.csv']));
assert(~isfile(fullfile(dest,'simulation_summary.csv')),'P5 already exists');
cfg=struct('smoke_trials',10,'pilot_trials',100,'starts',4,'extra_n2',4, ...
 'solver_seed',40912,'smoke_seed_base',410000,'pilot_seed_base',510000, ...
 'gain_threshold',.10,'min_area_fraction',.05,'min_separation_meV',4, ...
 'interpretation','predeclared engineering split flag, not statistically calibrated detection', ...
 'noise','hypothetical iid Gaussian; 1 percent max of scenario mean; no raw-frame bootstrap');
c20_v4_io('json',fullfile(dest,'P5_frozen_config.json'),cfg);
r=runtests('tests/C20V4SimulationTest.m'); writetable(table(r),fullfile(out,'tests','simulation_tests.csv')); assert(all([r.Passed]));
q=load(fullfile(out,'590_PL2_10w','centered_binned_spectra.mat'));
E=q.E(q.E>=300&q.E<=1800); members=q.bins(2).q_members;
scenarios={'static_single','q_mixed_single','overlapping_double','background_mismatch'};
rows={}; summaries={};
for phase=1:2
 count=cfg.smoke_trials; base=cfg.smoke_seed_base; if phase==2, count=cfg.pilot_trials; base=cfg.pilot_seed_base; end
 for si=1:4
  checkpoint=fullfile(dest,sprintf('simulation_phase%d_%s.mat',phase,scenarios{si}));
  trials=cell(1,count); localrows={};
  if isfile(checkpoint), previous=load(checkpoint); assert(isequal(previous.cfg,cfg)); trials=previous.trials; end
  for trial=1:count
   seed=base+si*1000+trial;
   if isempty(trials{trial})
    sim=c20_v4_simulate(E,members,scenarios{si},seed);
    fits=qe_compare_component_models(E,sim.observed,n_starts=cfg.starts,start_policy='independent', ...
     extra_starts=cfg.extra_n2,seed=cfg.solver_seed);
   else
    sim=trials{trial}.simulation; fits=trials{trial}.fits; assert(sim.seed==seed);
   end
   success=all([fits.success]); gain=NaN; split=false; boundary=false; collapse=false; errorE=NaN;
   if success
    gain=(fits(1).normalized_sse-fits(2).normalized_sse)/max(fits(1).normalized_sse,eps);
    f=fits(2); areas=trapz(E,f.components); fraction=min(areas)/max(sum(areas),eps);
    separation=diff(f.parameters(:,1)); collapse=separation<cfg.min_separation_meV;
    split=gain>cfg.gain_threshold&&fraction>cfg.min_area_fraction&&~collapse;
    boundary=any(f.candidates(f.selected_start).boundary);
    if si==3, errorE=max(abs(f.parameters(:,1)-sim.truth_parameters(:,1))); end
   end
   for n=1:2
    for k=1:numel(fits(n).candidates), fits(n).candidates(k).jacobian=[]; end
   end
   trials{trial}=struct('simulation',sim,'fits',fits);
   row={phase,string(scenarios{si}),trial,seed,success,gain,split,boundary,collapse,errorE};
   rows(end+1,:)=row; localrows(end+1,:)=row; %#ok<AGROW>
   if mod(trial,10)==0, disp(sprintf('P5 phase%d %s %d/%d',phase,scenarios{si},trial,count)); end
  end
  if ~isfile(checkpoint), save(checkpoint,'trials','cfg','-v7'); end
  L=cell2table(localrows); k=nnz(L{:,7}); phat=k/count; z=1.95996398454;
  center=(phat+z^2/(2*count))/(1+z^2/count); radius=z*sqrt(phat*(1-phat)/count+z^2/(4*count^2))/(1+z^2/count);
  summaries(end+1,:)={phase,string(scenarios{si}),count,nnz(~L{:,5}),k,phat,center-radius,center+radius,nnz(L{:,8}),nnz(L{:,9})}; %#ok<AGROW>
  if phase==1, assert(all(L{:,5}),'Smoke optimizer failures: pilot expansion stopped'); end
 end
end
writetable(cell2table(rows,'VariableNames',{'phase','scenario','trial','seed','success','relative_gain','engineering_split','n2_boundary','collapse','max_E0_error_meV'}),fullfile(dest,'simulation_trials.csv'));
writetable(cell2table(summaries,'VariableNames',{'phase','scenario','trials','failures','engineering_splits','fraction','Wilson_low','Wilson_high','boundary','collapse'}),fullfile(dest,'simulation_summary.csv'));
c20_v4_io('text',fullfile(dest,'INTERPRETATION.md'), ...
 'Simulations are hypothetical engineering tests, not measured C20 false-positive rates. The fixed split rule is uncalibrated. Rejected static H0 cannot rule out momentum or time mixing. No experimental mean or bootstrap residual was used in H0 generation. Full trial arrays and candidates are saved; Jacobian arrays omitted for packet size, singular values retained.');
end
