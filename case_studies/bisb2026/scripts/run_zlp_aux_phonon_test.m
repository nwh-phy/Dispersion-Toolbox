function out = run_zlp_aux_phonon_test()
% Does the 50-150 meV loss band carry a MoS2 optical-phonon component?
% A1 spectra from the v7 mirror-q run, DL n=2 B1 peaks, three aux layouts:
%   aux1   one free aux peak 30-300 meV (pilot setting)
%   aux2c  aux held to 44-56 meV (MoS2 E', A1', 2LA) + free aux 60-300 meV
%   aux2f  free aux 30-80 meV + free aux 80-300 meV (does the low one find ~50?)
% Each at the default ZLP core half-width (2.5 FWHM) and at 40 meV, which
% gives 45-65 meV full data weight. Parent outputs stay read-only.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p7=fullfile(root,'paper_results','b1_mirror_q_v7','20261007T074250Z');
cfg=struct('labels',{{'R1','R2','R3','R1m','R2m','R3m'}},'targets',[-.0025 .0075 .0125 .0025 -.0075 -.0125], ...
 'seed_offset',[1 2 3 1 2 3],'window',[300 1800],'model','lorentz','n_peaks',2,'J_starts',12, ...
 'seed',20260912,'core_halfwidths',[NaN 40]);
variants={'aux1',[30 300];'aux2c',[44 56;60 300];'aux2f',[30 80;80 300]};
id=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmss''Z'''));
out=fullfile(root,'paper_results','zlp_aux_phonon_test',[id '_A1']); assert(~isfolder(out));
mkdir(fullfile(out,'figures')); diary(fullfile(out,'execution.log')); dc=onCleanup(@()diary('off'));
disp(['RUN_DIR=' out]);
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root)); cfg.git_head=strtrim(head);
cfg.variants=variants;

v7=load(fullfile(p7,'mirror_fits.mat')); keys=cellfun(@(a)a.key,v7.A,'UniformOutput',false);
rows={}; R=struct([]);
for r=1:numel(cfg.labels)
 d7=v7.A{strcmp(keys,[cfg.labels{r} '_' cfg.model])}; E=d7.E1(:); Y=d7.Y1(:);
 seed=cfg.seed+cfg.seed_offset(r);
 for ch=cfg.core_halfwidths
  F=struct();
  for v=1:size(variants,1)
   fj=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
    'peak_model',cfg.model,'n_peaks',cfg.n_peaks,'n_starts',cfg.J_starts,'seed',seed, ...
    'aux_windows',variants{v,2},'core_halfwidth',ch);
   rows(end+1,:)=summary_row(cfg.labels{r},cfg.targets(r),ch,variants{v,1},fj); %#ok<AGROW>
   F.(variants{v,1})=fj;
  end
  R(end+1).label=cfg.labels{r}; R(end).q=cfg.targets(r); R(end).core_setting=ch; %#ok<AGROW>
  R(end).E=E; R(end).Y=Y; R(end).fits=F;
  plot_target(fullfile(out,'figures'),R(end));
 end
end
T=cell2table(rows,'VariableNames',{'region','q_Ainv','core_setting','core_halfwidth','variant','status', ...
 'cost','chi2_signal','chi2_low','chi2_gain','E0_1','E0_2', ...
 'a1_E0','a1_W','a1_height','a1_boundary','a2_E0','a2_W','a2_height','a2_boundary'});
writetable(T,fullfile(out,'aux_phonon_test.csv')); save(fullfile(out,'aux_phonon_test.mat'),'cfg','R','T','-v7.3');
disp(T);
end

function row=summary_row(label,q,ch,variant,f)
% chi2_low: loss side 30-300 meV from the ZLP centre, plain 1/sigma^2.
a=nan(2,4); c=NaN; cs=NaN; cl=NaN; cg=NaN; e0=[NaN NaN];
if f.success
 c=f.cost; cs=f.chi2_red_signal; cg=f.chi2_red_gain; e0=f.parameters(:,1).';
 x=f.energy_meV-f.zlp_parameters.center_meV; lo=x>=30&x<=300&isfinite(f.residual);
 cl=mean((f.residual(lo)./f.sigma(lo)).^2);
 for k=1:size(f.aux_parameters,1)
  a(k,:)=[f.aux_parameters(k,1:2) max(f.aux_peaks(:,k)) any(f.aux_boundary(k,:))];
 end
end
setting='default'; if isfinite(ch), setting=sprintf('%g',ch); end
row={label,q,setting,f.effective_options.core_halfwidth,variant,f.numerical_status,c,cs,cl,cg, ...
 e0(1),e0(2),a(1,1),a(1,2),a(1,3),a(1,4),a(2,1),a(2,2),a(2,3),a(2,4)};
end

function plot_target(dir,R)
f=figure('Visible','off','Position',[100 100 1650 480]); tiledlayout(1,3);
v=fieldnames(R.fits); hw=NaN;
for i=1:numel(v)
 j=R.fits.(v{i}); nexttile; semilogy(R.E,R.Y,'k.','MarkerSize',5,'DisplayName','data'); hold on;
 ttl=[v{i} ': failed'];
 if j.success
  hw=j.effective_options.core_halfwidth;
  semilogy(j.energy_meV,j.zlp,'b-','DisplayName','ZLP');
  semilogy(j.energy_meV,sum(j.peaks,2),'m:','DisplayName','B1 (2 DL)');
  for k=1:size(j.aux_peaks,2)
   semilogy(j.energy_meV,j.aux_peaks(:,k),'--','LineWidth',1.3,'DisplayName',sprintf('aux %d',k));
  end
  semilogy(j.energy_meV,j.prediction,'r-','DisplayName','total');
  x=j.energy_meV-j.zlp_parameters.center_meV; lo=x>=30&x<=300;
  ttl=sprintf('%s: aux E0 %s meV, \\chi^2_{30-300} %.2f',v{i},mat2str(round(j.aux_parameters(:,1).')), ...
   mean((j.residual(lo)./j.sigma(lo)).^2,'omitnan'));
 end
 xline(hw,'k:','HandleVisibility','off');
 xlim([-180 400]); ylim([max(min(R.Y(R.Y>0)),1) max(R.Y)*1.5]); grid on; title(ttl); legend('Location','northeast');
end
sgtitle(sprintf('%s q=%+.4f, DL n=2, core half-width %.0f meV',R.label,R.q,hw));
exportgraphics(f,fullfile(dir,sprintf('%s_core%.0f.png',R.label,hw)),'Resolution',120); close(f);
end
