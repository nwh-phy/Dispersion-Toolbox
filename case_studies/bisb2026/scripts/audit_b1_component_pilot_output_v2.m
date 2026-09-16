function audit_b1_component_pilot_output_v2(out)
% Independent post-run checks on saved arrays, not just solver return codes.
out=char(out); p=fullfile(out,'590_PL2_10w');
saved=load(fullfile(p,'fit_details.mat')); a=load(fullfile(p,'binned_spectra.mat'));
rows=table();
for d=saved.details
 d=d{1};
 for f=d.fits
  if ~f.success, continue; end
  model=peak_models(f.peak_model); curves=zeros(size(f.components));
  for j=1:f.n_components
   pp=f.parameters(j,:); curves(:,j)=model.model_fn(pp(1),pp(2),pp(3),f.energy_meV);
  end
  curve_error=max(abs(curves-f.components),[],'all');
  total_error=max(abs(f.background+sum(curves,2)-f.prediction));
  residual_error=max(abs(f.observed-f.prediction-f.residual));
  sse_error=abs(sum(f.residual.^2)-f.sse);
  assert(curve_error<1e-8*max(1,max(abs(f.observed))));
  assert(total_error<1e-8*max(1,max(abs(f.observed)))&&residual_error<1e-10&&sse_error<1e-8);
  rows=[rows;table(string(d.key),f.n_components,curve_error,total_error,residual_error,sse_error, ...
   'VariableNames',{'key','n','native_curve_error','total_error','residual_error','sse_error'})]; %#ok<AGROW>
 end
end
writetable(rows,fullfile(out,'saved_output_contract_checks.csv'));
for bi=1:numel(a.allbins)
 b=a.allbins{bi};
 for j=1:numel(b.units)
  u=b.units(j); assert(all(sign(str2double(split(string(u.source_q_Ainv),',')))==sign(u.q_Ainv)));
  assert(all(diff(u.source_channel)==1));
  assert(max(abs(b.mean(:,j)*u.source_q_count-b.sum(:,j)))<1e-10*max(1,max(abs(b.sum(:,j)))));
 end
end
% Re-run one representative pair with exactly the frozen start schedule.
d=saved.details{1}; original=d.fits; cfg=jsondecode(fileread(fullfile(out,'config_resolved.json')));
again=qe_compare_component_models(original(1).energy_meV,original(1).observed, ...
 energy_window=[min(original(1).energy_meV) max(original(1).energy_meV)], ...
 peak_model=original(1).peak_model,n_starts=cfg.n_starts,seed=cfg.seed);
for n=1:2
 assert(max(abs(again(n).prediction-original(n).prediction))<1e-8*max(original(n).observed));
end
fid=fopen(fullfile(out,'reproducibility_check.txt'),'w');
fprintf(fid,'Saved native-parameter curve reconstruction, total/residual/SSE contracts and one representative 24-start pair rerun passed.\n'); fclose(fid);
end
