function b1_fullq_stage0c(out, labels)
% B1 full-q task, stage 0c: recheck of the low-q degeneracy with the
% prefactor held below cfg.prefactor_floor and not applied to aux peaks.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
p7=fullfile(root,'paper_results','b1_mirror_q_v7','20261007T074250Z');
cfg=b1_fullq_config(); all_labels={'R1','R2','R3','R1m','R2m','R3m'}; seed_offset=[1 2 3 1 2 3];
l=load(fullfile(p2,'590_PL2_10w','L1_minimal.mat')); qe=l.d.qe; sp=4.12e-4;
v7=load(fullfile(p7,'mirror_fits.mat')); keys=cellfun(@(a)a.key,v7.A,'UniformOutput',false);
if ~isfolder(out), mkdir(out); end
V={'none','',NaN,0,true;'2d_h20_floor','2d',0.002,cfg.prefactor_floor,false;'2d_h5_floor','2d',0.0005,cfg.prefactor_floor,false; ...
   '2d_h80_floor','2d',0.008,cfg.prefactor_floor,false;'3d_h20_floor','3d',0.002,cfg.prefactor_floor,false;'2d_h20_nofloor','2d',0.002,0,true};
for li=1:numel(labels)
 r=find(strcmp(all_labels,labels{li})); d7=v7.A{strcmp(keys,[labels{li} '_lorentz'])};
 E=d7.E1(:); Y=d7.Y1(:); mq=qe.q_Ainv(d7.n3_channels); rows={};
 for v=1:size(V,1)
  pre=[]; if ~isempty(V{v,2}), pre=qe_kinematic_prefactor(mq,form=V{v,2},h_perp=V{v,3},sigma_probe=sp); end
  f=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf 1800],'signal_window',[300 1800],'peak_model','lorentz','n_peaks',2, ...
   'n_starts',12,'seed',20260912+seed_offset(r),'aux_windows',cfg.aux{1,2},'prefactor',pre, ...
   'prefactor_floor_meV',V{v,4},'prefactor_on_aux',V{v,5});
  x=f.energy_meV-f.zlp_parameters.center_meV; lo=x>=30&x<=300; Ep=f.energy_meV>0;
  a=trapz(f.energy_meV(Ep),f.loss_function(Ep,:));
  rows(end+1,:)={labels{li},V{v,1},f.cost,f.chi2_red_signal,mean((f.residual(lo)./f.sigma(lo)).^2),f.chi2_red_gain, ...
   f.parameters(1,1),f.parameters(1,2),f.parameters(2,1),f.parameters(2,2),a(1)/sum(a),f.numerical_status}; %#ok<AGROW>
  fprintf('%s %-15s cost %8.1f chi2 sig %.3f low %6.2f gain %.2f | E0 %.0f/%.0f W %.0f/%.0f w %.3f %s\n',rows{end,1:2}, ...
   rows{end,3},rows{end,4},rows{end,5},rows{end,6},rows{end,7},rows{end,9},rows{end,8},rows{end,10},rows{end,11},rows{end,12});
 end
 writetable(cell2table(rows,'VariableNames',{'region','variant','cost','chi2_signal','chi2_low','chi2_gain','E0_1','W_1', ...
  'E0_2','W_2','weight_low','status'}),fullfile(out,sprintf('stage0c_%s.csv',labels{li})));
end
end
