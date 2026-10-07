function b1_fullq_hscan(out, labels)
% B1 full-q task, stage 0b: chi2 scan of the unknown perpendicular acceptance
% half-width h in the kinematic prefactor (2D and 3D forms), DL n=2, main aux,
% on v7 A1 spectra. Writes hscan_<label>.csv into out.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
p7=fullfile(root,'paper_results','b1_mirror_q_v7','20261007T074250Z');
pp=fullfile(root,'paper_results','b1_physics_preview','20260912T064410245Z_5756ece8_93247099');
all_labels={'R1','R2','R3','R1m','R2m','R3m'}; seed_offset=[1 2 3 1 2 3];
hs=[0 0.001 0.002 0.003 0.005 0.008 0.012 0.02];
l=load(fullfile(p2,'590_PL2_10w','L1_minimal.mat')); qe=l.d.qe;
m=load(fullfile(pp,'590_PL2_10w','A1_full_q_map.mat')); map=m.map;
el=sum(map.A1(abs(map.E)<=12,:),1); el(map.invalid_q)=0; [~,i0]=max(el); k=-2:2; w=el(i0+k); w=w/sum(w);
sp=sqrt(max(sum(w.*(k-sum(w.*k)).^2)-1/12,0))*median(diff(map.q));
v7=load(fullfile(p7,'mirror_fits.mat')); keys=cellfun(@(a)a.key,v7.A,'UniformOutput',false);
if ~isfolder(out), mkdir(out); end
for li=1:numel(labels)
 r=find(strcmp(all_labels,labels{li})); d7=v7.A{strcmp(keys,[labels{li} '_lorentz'])};
 E=d7.E1(:); Y=d7.Y1(:); mq=qe.q_Ainv(d7.n3_channels); rows={};
 V=[{'none',NaN};[repmat({'2d'},numel(hs),1) num2cell(hs(:))];[repmat({'3d'},numel(hs),1) num2cell(hs(:))]];
 for v=1:size(V,1)
  pre=[]; if ~strcmp(V{v,1},'none'), pre=qe_kinematic_prefactor(mq,form=V{v,1},h_perp=V{v,2},sigma_probe=sp); end
  f=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf 1800],'signal_window',[300 1800],'peak_model','lorentz', ...
   'n_peaks',2,'n_starts',12,'seed',20260912+seed_offset(r),'aux_windows',[30 80;80 300],'prefactor',pre);
  x=f.energy_meV-f.zlp_parameters.center_meV; lo=x>=30&x<=300; Ep=f.energy_meV>0;
  a=trapz(f.energy_meV(Ep),f.loss_function(Ep,:));
  rows(end+1,:)={labels{li},mean(mq),V{v,1},V{v,2},f.cost,f.chi2_red_signal,mean((f.residual(lo)./f.sigma(lo)).^2), ...
   f.chi2_red_gain,f.parameters(1,1),f.parameters(1,2),f.parameters(2,1),f.parameters(2,2),a(1)/sum(a),f.numerical_status}; %#ok<AGROW>
  fprintf('%s %s h=%g: cost %.1f chi2 sig %.3f low %.2f | E0 %.0f/%.0f W %.0f/%.0f wlow %.3f\n',labels{li},V{v,1},V{v,2}, ...
   rows{end,5},rows{end,6},rows{end,7},rows{end,9},rows{end,11},rows{end,10},rows{end,12},rows{end,13});
 end
 writetable(cell2table(rows,'VariableNames',{'region','q_mean','form','h','cost','chi2_signal','chi2_low','chi2_gain', ...
  'E0_1','W_1','E0_2','W_2','weight_low','status'}),fullfile(out,sprintf('hscan_%s.csv',labels{li})));
end
end
