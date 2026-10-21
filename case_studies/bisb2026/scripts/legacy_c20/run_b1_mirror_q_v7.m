function out = run_b1_mirror_q_v7()
% C20/B1 +/-q mirror check. Re-runs the v5 A1 N=3 two-component fits and the
% five-member M1/M2 comparison at the three v5 targets and their mirror points,
% with identical settings, so each +q/-q pair is compared on equal footing.
% Parent outputs stay read-only.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
p5=fullfile(root,'paper_results','b1_components_v5','20260912T055651568Z_c1dfab1b_fcd805b3');
raw_path=fullfile(root,'20260120 BiSb','590 PL2 10w 0.004 10sx300','Sequence EELS Image 274.npy');
% Mirror targets reuse the seed of their v5 partner region.
cfg=struct('targets',[-.0025 .0075 .0125 .0025 -.0075 -.0125],'labels',{{'R1','R2','R3','R1m','R2m','R3m'}}, ...
 'seed_offset',[1 2 3 1 2 3],'window',[300 1800],'reference_window',[300 1800], ...
 'A_starts',24,'A_extra_n2',12,'member_starts',12,'seed',20260912, ...
 'models',{{'lorentz_symmetric','lorentz'}},'zlp_window',[-100 100],'reference_abs_q',.001);
id=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmss''Z'''));
out=fullfile(root,'paper_results','b1_mirror_q_v7',id); assert(~isfolder(out));
mkdir(fullfile(out,'figures')); diary(fullfile(out,'execution.log')); dc=onCleanup(@()diary('off'));
disp(['RUN_DIR=' out]);
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root));
c20_v4_io('json',fullfile(out,'config.json'),struct('cfg',cfg,'git_head',strtrim(head),'parent_v2',p2,'parent_v5',p5, ...
 'matlab',version,'arch',computer('arch')));

l=load(fullfile(p2,'590_PL2_10w','L1_minimal.mat')); qe=l.d.qe; E=qe.energy_meV;
n3=qe_centered_bins(qe,cfg.targets,3); n5=qe_centered_bins(qe,cfg.targets,5);
assert(all([n3.valid])&&all([n5.valid]),'Invalid centered bins');

raw=read_npy(raw_path); T=size(raw,1);
assert(size(raw,2)==numel(qe.q_Ainv)&&size(raw,3)==numel(E),'Raw axes mismatch');
invalid=false(1,size(raw,2));
for q=1:size(raw,2)
 x=raw(:,q,:); invalid(q)=any(~isfinite(x(:))|double(x(:))>=double(intmax('uint32')));
end
assert(isequal(find(invalid),find(any(isnan(qe.intensity),1))),'Invalid mask differs from parent');
members=unique([[n3.source_channel] [n5.source_channel]]); assert(~any(invalid(members)));
X=permute(double(raw(:,members,:)),[3 1 2]);                       % E x T x member
refq=find(abs(qe.q_Ainv)<=cfg.reference_abs_q & ~invalid);
reference=squeeze(sum(double(raw(:,refq,:)),2)).'; clear raw
al=qe_zlp_integer_align(E,reference,X,cfg.zlp_window);
fprintf('A1 alignment: %d/%d valid frames, offsets %d..%d px\n',nnz(al.valid),T, ...
 min(al.measured_offset_pixels(al.valid)),max(al.measured_offset_pixels(al.valid)));

A={}; M={}; rows={};
for r=1:numel(cfg.targets)
 b3=n3(r); b5=n5(r); [~,i3]=ismember(b3.source_channel,members); [~,i5]=ismember(b5.source_channel,members);
 Y1=sum(mean(al.aligned(:,:,i3),3),2);                               % A1, sequence-summed N=3 mean
 Xm=squeeze(sum(X(al.support,:,i5),2));                              % A0 five members, as in v5
 seed=cfg.seed+cfg.seed_offset(r);
 for mi=1:numel(cfg.models)
  model=cfg.models{mi};
  fa=qe_compare_component_models(al.E,Y1,energy_window=cfg.window,peak_model=model, ...
   n_starts=cfg.A_starts,extra_starts=cfg.A_extra_n2,start_policy='independent',seed=seed);
  f0=qe_compare_component_models(E,b3.mean,energy_window=cfg.window,peak_model=model, ...
   n_starts=cfg.A_starts,extra_starts=cfg.A_extra_n2,start_policy='independent',seed=seed);
  fm=qe_fit_member_models(al.E,qe.q_Ainv(b5.source_channel),Xm,energy_window=cfg.window, ...
   peak_model=model,n_starts=cfg.member_starts,seed=seed,initial_fits=f0);
  key=sprintf('%s_%s',cfg.labels{r},model);
  d=struct('key',key,'target',cfg.targets(r),'n3_channels',b3.source_channel,'n5_channels',b5.source_channel, ...
   'A1',fa,'A0_init',f0,'members',fm,'Y1',Y1,'E1',al.E);
  A{end+1}=d; %#ok<AGROW>
  c20_v4_plot_fits(struct('key',['A1_' key],'unit',struct('q_Ainv',b3.q_Ainv,'source_channel',b3.source_channel),'fits',fa), ...
   fullfile(out,'figures',['A1_' key '.png']));
  f2=fa(2); apex=nan(1,2); area=nan(1,2);
  if f2.success
   for j=1:2, [~,k]=max(f2.components(:,j)); apex(j)=f2.energy_meV(k); end
   s=qe_area_scopes(f2.energy_meV,f2.components,cfg.reference_window); area=s.area_reference_window;
  end
  red=100*(1-fm(2).objective/fm(1).objective);
  rows(end+1,:)={string(cfg.labels{r}),cfg.targets(r),string(model),string(mat2str(b3.source_channel)), ...
   apex(1),apex(2),area(1)/sum(area),f2.parameters(1,2),f2.parameters(2,2),string(f2.numerical_status), ...
   100*(1-fa(2).normalized_sse/fa(1).normalized_sse),fm(1).objective,fm(2).objective,red, ...
   fm(2).slope_meV_A(1),fm(2).slope_meV_A(2)}; %#ok<AGROW>
  fprintf('%-26s apex %4.0f / %4.0f meV  low frac %.3f  A1 n2 drop %.1f%%  member M2 drop %.1f%%\n', ...
   key,apex(1),apex(2),area(1)/sum(area),rows{end,11},red);
 end
end
T=cell2table(rows,'VariableNames',{'region','target_q_Ainv','peak_model','n3_channels', ...
 'apex_low_meV','apex_high_meV','low_area_fraction_W300_1800','width_low_meV','width_high_meV','A1_n2_status', ...
 'A1_n2_vs_n1_objective_drop_pct','member_M1_objective','member_M2_objective','member_M2_vs_M1_drop_pct', ...
 'member_M2_slope_low_meV_per_Ainv','member_M2_slope_high_meV_per_Ainv'});
writetable(T,fullfile(out,'mirror_summary.csv'));
save(fullfile(out,'mirror_fits.mat'),'A','cfg','-v7');

% Reproduction check against v5 for the three original regions.
v5=load(fullfile(p5,'member_models','all_member_fits.mat')); chk={};
for k=1:numel(v5.allmembers)
 d5=v5.allmembers{k}; ix=find(strcmp(cellfun(@(a)a.key,A,'UniformOutput',false),d5.key),1);
 a5=load(fullfile(p5,'A0_A1',['A1_' d5.key '.mat']));
 chk(end+1,:)={string(d5.key),d5.fits(2).objective,A{ix}.members(2).objective, ...
  max(abs(a5.d.fits(2).parameters(:,1)-A{ix}.A1(2).parameters(:,1)))}; %#ok<AGROW>
end
C=cell2table(chk,'VariableNames',{'key','v5_member_M2_objective','v7_member_M2_objective','A1_n2_max_center_diff_meV'});
writetable(C,fullfile(out,'v5_reproduction_check.csv')); disp(C);

% +/-q overlay of the A1 spectra and two-component fits.
pairs=[1 4; 2 5; 3 6];
for mi=1:numel(cfg.models)
 f=figure('Visible','off','Position',[100 100 1500 480]); tiledlayout(1,3);
 for p=1:3
  nexttile; hold on;
  for s=1:2
   d=A{(pairs(p,s)-1)*numel(cfg.models)+mi}; z=d.A1(2); c=lines(2);
   plot(z.energy_meV,z.observed/max(z.observed),'.','Color',c(s,:),'MarkerSize',4,'DisplayName',sprintf('q=%+.4f obs',d.target));
   plot(z.energy_meV,z.components/max(z.observed),'-','Color',c(s,:),'LineWidth',1.2,'HandleVisibility','off');
  end
  xlim(cfg.window); xlabel('Energy loss (meV)'); ylabel('A1 counts / max'); legend('Location','northeast');
  title(sprintf('|q|=%.4f 1/A, %s',abs(cfg.targets(pairs(p,1))),cfg.models{mi}),'Interpreter','none');
 end
 exportgraphics(f,fullfile(out,'figures',sprintf('pm_q_overlay_%s.png',cfg.models{mi}))); close(f);
end
disp(['DONE=' out]);
end
