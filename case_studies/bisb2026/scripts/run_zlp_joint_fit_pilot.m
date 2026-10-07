function out = run_zlp_joint_fit_pilot()
% Pilot: ZEST-style joint ZLP fit vs the current power-law window fit on the
% v7 mirror-q targets (590, A0 N=3 means). Same spectra, same n=1/2
% components; only the background treatment differs. Joint variants:
%   full  all of [-180, 1800] meV in the data term
%   gap   loss side 30-300 meV dropped; ZLP tail set by gain side + core
%   aux   one detailed-balance DL peak in 30-300 meV for the low-energy losses
% Parent outputs stay read-only.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
cfg=struct('targets',[-.0025 .0075 .0125 .0025 -.0075 -.0125],'labels',{{'R1','R2','R3','R1m','R2m','R3m'}}, ...
 'seed_offset',[1 2 3 1 2 3],'window',[300 1800],'models',{{'lorentz_symmetric','lorentz'}}, ...
 'A_starts',24,'A_extra_n2',12,'J_starts',12,'seed',20260912,'low_window',[30 300]);
id=char(datetime('now','TimeZone','UTC','Format','yyyyMMdd''T''HHmmss''Z'''));
out=fullfile(root,'paper_results','zlp_joint_fit_pilot',id); assert(~isfolder(out));
mkdir(fullfile(out,'figures')); diary(fullfile(out,'execution.log')); dc=onCleanup(@()diary('off'));
disp(['RUN_DIR=' out]);
[~,head]=system(sprintf('git -C "%s" rev-parse HEAD',root)); cfg.git_head=strtrim(head);

l=load(fullfile(p2,'590_PL2_10w','L1_minimal.mat')); qe=l.d.qe; E=qe.energy_meV(:);
n3=qe_centered_bins(qe,cfg.targets,3); assert(all([n3.valid]),'Invalid centered bins');
variants={'full',{};'gap',{'exclude_windows',cfg.low_window};'aux',{'aux_windows',cfg.low_window}};

rows={}; R=struct([]);
for r=1:numel(cfg.targets)
 b=n3(r); Y=b.mean(:); seed=cfg.seed+cfg.seed_offset(r);
 for mi=1:numel(cfg.models)
  model=cfg.models{mi};
  fa=qe_compare_component_models(E,Y,energy_window=cfg.window,peak_model=model, ...
   n_starts=cfg.A_starts,extra_starts=cfg.A_extra_n2,start_policy='independent',seed=seed);
  for n=1:2
   rows(end+1,:)=summary_row(cfg.labels{r},b.q_Ainv,model,n,'powerlaw_window',fa(n),0); %#ok<AGROW>
   J=struct();
   for v=1:size(variants,1)
    t0=tic; fj=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
     'peak_model',model,'n_peaks',n,'n_starts',cfg.J_starts,'seed',seed,variants{v,2}{:});
    rows(end+1,:)=summary_row(cfg.labels{r},b.q_Ainv,model,n,['joint_' variants{v,1}],fj,toc(t0)); %#ok<AGROW>
    J.(variants{v,1})=fj;
   end
   R(end+1).label=cfg.labels{r}; R(end).q=b.q_Ainv; R(end).model=model; R(end).n=n; %#ok<AGROW>
   R(end).Y=Y; R(end).powerlaw=fa(n); R(end).joint=J;
   if n==2, plot_target(fullfile(out,'figures'),R(end),E,cfg.window); end
  end
 end
end
T=cell2table(rows,'VariableNames',{'region','q_Ainv','model','n','method','status','E0_1','W_1','E0_2','W_2', ...
 'bg_300','bg_400','bg_600','bg_fraction','chi2_red_signal','aux_E0','aux_W','seconds'});
writetable(T,fullfile(out,'comparison.csv')); save(fullfile(out,'pilot_results.mat'),'cfg','R','T','-v7.3');
disp(T);
end

function row=summary_row(label,q,model,n,method,f,secs)
p=nan(2,2); bg=nan(1,3); bf=NaN; chi=NaN; aux=[NaN NaN];
if f.success
 p(1:n,:)=f.parameters(1:n,1:2);
 Ef=f.energy_meV(:); bg=interp1(Ef,f.background(:),[300 400 600]);
 if isfield(f,'background_fraction_signal')
  bf=f.background_fraction_signal; chi=f.chi2_red_signal;
  if ~isempty(f.aux_parameters), aux=f.aux_parameters(1,1:2); end
 else
  in=Ef>=300&Ef<=1800; bf=sum(f.background(in))/sum(f.observed(in));
 end
end
row={label,q,model,n,method,f.numerical_status,p(1,1),p(1,2),p(2,1),p(2,2),bg(1),bg(2),bg(3),bf,chi,aux(1),aux(2),secs};
end

function plot_target(dir,R,E,win)
f=figure('Visible','off','Position',[100 100 1600 520]); tiledlayout(1,2);
c=struct('full',[0 .45 .74],'gap',[.47 .67 .19],'aux',[.85 .33 .1]); v=fieldnames(R.joint);
a=R.powerlaw; Y=R.Y;
nexttile; semilogy(E,Y,'k.','MarkerSize',4); hold on;
for i=1:numel(v)
 j=R.joint.(v{i}); if ~j.success, continue; end
 semilogy(j.energy_meV,j.zlp,'-','Color',c.(v{i}),'DisplayName',['ZLP ' v{i}]);
 if ~isempty(j.aux_peaks), semilogy(j.energy_meV,j.aux_peaks,'--','Color',c.(v{i}),'DisplayName','aux peak'); end
end
j=R.joint.gap; if j.success, semilogy(j.energy_meV,j.prefit_zlp,'k:','DisplayName','prefit ZLP'); end
if a.success && any(a.background>0), semilogy(a.energy_meV,a.background,'m--','DisplayName','power law'); end
xlim([E(1) win(2)]); ylim([max(min(Y(Y>0)),1) max(Y)*1.5]); grid on; legend('Location','northeast');
title(sprintf('%s q=%+.4f %s n=%d (log)',R.label,R.q,R.model,R.n),'Interpreter','none');
nexttile; m=E>=win(1)-100&E<=win(2); plot(E(m),Y(m),'k.','MarkerSize',5,'HandleVisibility','off'); hold on;
s={};
for i=1:numel(v)
 j=R.joint.(v{i}); if ~j.success, continue; end
 jm=j.energy_meV>=win(1)-100;
 plot(j.energy_meV(jm),j.prediction(jm),'-','Color',c.(v{i}),'DisplayName',['total ' v{i}]);
 plot(j.energy_meV(jm),j.background(jm)+j.peaks(jm,:),':','Color',c.(v{i}),'HandleVisibility','off');
 s{end+1}=sprintf('%s: %s',v{i},j.numerical_status); %#ok<AGROW>
end
if a.success, plot(a.energy_meV,a.background+a.components,'m:','HandleVisibility','off'); s{end+1}=['power law: ' a.numerical_status]; end
xline(win(1),'k-','HandleVisibility','off'); xlim([win(1)-100 win(2)]); grid on; legend('Location','northeast');
title(strjoin(s,' | '),'Interpreter','none');
exportgraphics(f,fullfile(dir,sprintf('%s_%s_n2.png',R.label,R.model)),'Resolution',120); close(f);
end
