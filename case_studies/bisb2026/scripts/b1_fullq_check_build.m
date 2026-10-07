function b1_fullq_check_build(run_root, session)
% Stage 1 check: frame-summed bins must equal the stored A1 full-q map
% (sum of the member columns) and, for 590 at +/-0.0025, 3x the v7 Y1.
root=bisb_find_project_root(fileparts(mfilename('fullpath')));
pp=fullfile(root,'paper_results','b1_physics_preview','20260912T064410245Z_5756ece8_93247099');
p7=fullfile(root,'paper_results','b1_mirror_q_v7','20261007T074250Z');
S=load(fullfile(run_root,'stage1',sprintf('frames_%s.mat',session)),'D'); D=S.D;
m=load(fullfile(pp,session,'A1_full_q_map.mat')); map=m.map;
[Ec,ia,ib]=intersect(round(D.E),round(map.E)); worst=0;
for b=1:numel(D.bins)
 ref=sum(map.A1(ib,D.bins(b).source_channel),2); x=D.Ysum(ia,b);
 worst=max(worst,max(abs(x-ref))/max(abs(ref)));
end
fprintf('%s: %d bins, common E %g..%g meV (%d pts); max |Ysum - A1 map| / max = %.2e\n',session,numel(D.bins),Ec(1),Ec(end),numel(Ec),worst);
if strcmp(session,'590_PL2_10w')
 v7=load(fullfile(p7,'mirror_fits.mat')); keys=cellfun(@(a)a.key,v7.A,'UniformOutput',false);
 for lab={'R1','R1m'}
  d7=v7.A{strcmp(keys,[lab{1} '_lorentz'])}; b=find(abs(D.q_center-d7.target)<1e-9);
  [~,i1,i2]=intersect(round(D.E),round(d7.E1)); r=max(abs(D.Ysum(i1,b)-3*d7.Y1(i2)))/max(3*d7.Y1(i2));
  fprintf('  %s q=%+.4f: max |Ysum - 3*v7 Y1| / max = %.2e\n',lab{1},d7.target,r);
 end
end
lossr=D.E>=300&D.E<=1800; gainr=D.E<=-60; zl=abs(D.E)<=30;
fprintf('  noise ratio var/mean (common mode removed): loss %.3f, gain %.3f, ZLP core %.3f; raw loss %.3f\n', ...
 median(D.noise_ratio(lossr,:),'all','omitnan'),median(D.noise_ratio(gainr,:),'all','omitnan'), ...
 median(D.noise_ratio(zl,:),'all','omitnan'),median(D.noise_ratio_raw(lossr,:),'all','omitnan'));
fprintf('  lag-1/2/3 correlation along E: loss %s, gain %s\n',mat2str(round(median(D.noise_rho_loss,2,'omitnan')',3)),mat2str(round(median(D.noise_rho_gain,2,'omitnan')',3)));
fprintf('  frame intensity drift (sd of frame total / mean): %.3f; ZLP offsets %d..%d px; sigma_probe %.2e\n', ...
 median(std(D.frame_drift,0,1)),min(D.offsets),max(D.offsets),D.sigma_probe);
end
