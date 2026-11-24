function out = bisb_slot_analyze(cfg)
%BISB_SLOT_ANALYZE DM4 spatial scan -> local elastic alignment -> slot spectra.
% Spatial scan rows, not independent time frames, are bootstrap units.
% No q folding, transport, kinematic division or assumed diffraction indexing.
if nargin<1, cfg=bisb_slot_config(); end
root=bisb_find_project_root(); addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'third_party','read_dm_file'));
d=bisb_datasets(cfg.dataset);
out=fullfile(root,'paper_results','bi_slot_20261117',char(datetime('now','Format','yyyyMMdd_HHmmss')));
assert(~isfolder(out),'Output already exists'); mkdir(out);
[m,de,units,meta]=ReadDMFile(d.path,fullfile(out,'dm4_metadata.txt'),1);
assert(strcmp(units,'eV') && isequal(size(m),[1028 512 19 19]));
assert(all(isfinite(m(:))) && all(m(:)>=0),'Invalid raw samples');
Eraw=((0:size(m,1)-1)'-double(meta{8}))*de*1000;
q=cfg.channels; raw=reshape(m,size(m,1),size(m,2),[]); clear m
raw_map=sum(double(raw),3); [~,p0]=max(sum(raw_map,2)); origin=Eraw(p0);
% One ZLP shift per spatial position, measured on the strong direct beam and
% applied to every q channel. Independent weak-channel alignment locks onto noise.
elastic0=sum(raw_map(abs(Eraw-origin)<=100,:),1); [~,jd]=max(elastic0(q)); direct=q(jd);
refq=direct+(-2:2); ref=squeeze(sum(double(raw(:,refq,:)),2));
X=permute(raw(:,q,:),[1 3 2]);
al=qe_zlp_integer_align(Eraw-origin,ref,X,[-100 100]); clear X ref
assert(all(al.valid),'Unreliable direct-beam references require review');
E=Eraw(al.support)-Eraw(al.reference_pixel);
F=permute(al.aligned,[1 3 2]); al=rmfield(al,'aligned');
M=sum(F,3); offsets=al.measured_offset_pixels*de*1000;
% Per-channel sums over each scan row (19 consecutive positions) for bisb_slot_geometry.
Fr=squeeze(sum(reshape(F,size(F,1),size(F,2),19,19),3));
% Compare sum-then-align with spatial-position alignment: exact spectrum at each common pixel.
static=zeros(size(M));
for j=1:numel(q)
 y=raw_map(:,q(j)); ii=find(abs(Eraw-origin)<220); [~,p]=max(y(ii)); p=ii(p);
 ix=(1:numel(Eraw))'-p; static(:,j)=interp1(ix*de*1000,y,E,'linear',NaN);
end
% Common energy support has no circular wrapping; verify mapping on one spectrum.
j=17; t=31; ix=al.support+al.measured_offset_pixels(t);
assert(isequal(F(:,j,t),double(raw(ix,q(j),t))));
clear raw
nb=ceil(numel(q)/cfg.bin_width); S=zeros(numel(E),nb,361); Sb=zeros(numel(E),nb); qb=zeros(1,nb);
for b=1:nb
 jj=(b-1)*cfg.bin_width+1:min(b*cfg.bin_width,numel(q));
 S(:,b,:)=sum(F(:,jj,:),2); Sb(:,b)=sum(static(:,jj),2); qb(b)=mean(q(jj));
end
clear F
B=sum(S,3); elastic=sum(M(abs(E)<=40,:),1); zwidth=zeros(size(q));
for j=1:numel(q), [~,zwidth(j)]=qe_zlp_center(M(:,j),E); end
% Peak energies from each background choice; confidence intervals from spatial rows.
rng(cfg.seed); nw=size(cfg.background_windows,1); peaks=nan(nb,nw); exponent=peaks;
R=zeros(size(B)); boot=nan(nb,cfg.bootstrap); rawboot=boot; static_peak=nan(nb,1); raw_peak=static_peak;
for b=1:nb
 for w=1:nw
  [bg,r]=qe_powerlaw_tail(E,B(:,b),cfg.background_windows(w,:)); exponent(b,w)=r;
  if r>0, peaks(b,w)=qe_peak_apex(E,B(:,b)-bg,cfg.peak_window,smooth=15,half_width=60); end
  if w==2, R(:,b)=B(:,b)-bg; end
 end
 raw_peak(b)=qe_peak_apex(E,B(:,b),cfg.peak_window,smooth=15,half_width=60);
 [bg,r]=qe_powerlaw_tail(E,Sb(:,b),cfg.background_windows(2,:));
 if r>0, static_peak(b)=qe_peak_apex(E,Sb(:,b)-bg,cfg.peak_window,smooth=15,half_width=60); end
 rows=squeeze(sum(reshape(S(:,b,:),numel(E),19,19),2));
 for k=1:cfg.bootstrap
  y=sum(rows(:,randi(19,1,19)),2); [bg,r]=qe_powerlaw_tail(E,y,cfg.background_windows(2,:));
  if r>0, boot(b,k)=qe_peak_apex(E,y-bg,cfg.peak_window,smooth=15,half_width=60); end
  rawboot(b,k)=qe_peak_apex(E,y,cfg.peak_window,smooth=15,half_width=60);
 end
end
ci=prctile(boot,[2.5 97.5],2); success=sum(isfinite(boot),2)/cfg.bootstrap;
rawci=prctile(rawboot,[2.5 97.5],2); rawsuccess=sum(isfinite(rawboot),2)/cfg.bootstrap;
T=table(qb',raw_peak,rawci(:,1),rawci(:,2),rawsuccess,peaks(:,1),peaks(:,2),peaks(:,3),ci(:,1),ci(:,2),success,static_peak,exponent(:,2), ...
 'VariableNames',{'channel','raw_apex_meV','raw_boot_lo','raw_boot_hi','raw_success','background1_meV','background2_meV','background3_meV','row_boot_lo','row_boot_hi','bootstrap_success','sum_then_align_meV','tail_exponent'});
writetable(T,fullfile(out,'peak_sensitivity.csv'));
save(fullfile(out,'analysis.mat'),'cfg','d','E','q','M','Fr','qb','B','S','R','static','raw_map','Eraw','offsets','zwidth','elastic','T','al','meta','rawboot','boot','-v7.3');
qe_plot_heatmap(E,qb,B,fullfile(out,'01_slot_dispersion'),q_window=[min(qb) max(qb)],e_window=[200 3000],norm_window=[400 3000],q_label='Slot detector channel (uncalibrated)',title='Bi: 8-channel bins; local ZLP alignment');
qe_plot_heatmap(E,qb,B,fullfile(out,'02_low_energy'),q_window=[min(qb) max(qb)],e_window=[200 1500],norm_window=[400 1500],q_label='Slot detector channel (uncalibrated)',title='Bi: measured spectra, no background subtraction');
f=figure('Visible','off','Color','w','Position',[50 50 1200 850]); tiledlayout(2,2);
nexttile; imagesc(1:512,Eraw/1000,log10(max(raw_map,1))); axis xy; ylim([-.05 .45]); title('Raw elastic line'); xlabel('Detector channel'); ylabel('Raw energy (eV)'); colorbar;
nexttile; plot(1:361,offsets,'.-'); xlabel('Scan position in acquisition order'); ylabel('ZLP shift (meV)'); grid on; title('Direct-beam shift per spatial position');
nexttile; semilogy(q,elastic); xlabel('Detector channel'); ylabel('Counts within +/-40 meV'); grid on;
nexttile; plot(q,zwidth); xlabel('Detector channel'); ylabel('Aligned ZLP FWHM (meV)'); grid on;
exportgraphics(f,fullfile(out,'03_elastic_diagnostics.png'),'Resolution',180); close(f);
centers=[72 96 136 184 224 232 240 280 328 376];
f=figure('Visible','off','Color','w','Position',[50 50 1400 850]); tiledlayout(2,5);
for c=centers
 [~,b]=min(abs(qb-c)); nexttile; y=B(:,b); bg=y-R(:,b); scale=max(y(E>=400&E<=1200));
 plot(E/1000,y/scale,'k',E/1000,bg/scale,'--',E/1000,R(:,b)/scale,'b'); xlim([.15 1.6]); ylim([-.3 2]); grid on;
 title(sprintf('Channels %.0f-%.0f',qb(b)-3.5,min(qb(b)+3.5,max(q)))); xlabel('Loss (eV)');
end
legend('Measured','Power-law tail','Residual','Location','best');
exportgraphics(f,fullfile(out,'04_spectra.png'),'Resolution',180); close(f);
f=figure('Visible','off','Color','w','Position',[50 50 1100 650]);
errorbar(qb,raw_peak,raw_peak-rawci(:,1),rawci(:,2)-raw_peak,'o'); hold on; plot(qb,peaks,'.');
xlabel('Slot detector channel'); ylabel('Peak apex (meV)'); ylim([350 1250]); grid on;
labels={'Measured; 95% spatial-row bootstrap'};
for w=1:nw, labels{end+1}=sprintf('%g-%g meV tail',cfg.background_windows(w,:)); end
legend(labels,'Location','best');
exportgraphics(f,fullfile(out,'05_peak_sensitivity.png'),'Resolution',180); close(f);
fid=fopen(fullfile(out,'summary.txt'),'w'); clean=onCleanup(@()fclose(fid));
fprintf(fid,'Dataset %s\nVersion %s\nEnergy step %.9g meV\n%d/%d local ZLP references valid\nLocal shifts range %.1f to %.1f meV\nMedian aligned ZLP FWHM %.2f meV (5-95%% %.2f-%.2f)\n', ...
 d.id,cfg.version,de*1000,sum(al.valid),numel(al.valid),min(offsets(:)),max(offsets(:)),median(zwidth),prctile(zwidth,5),prctile(zwidth,95));
fprintf(fid,'No reciprocal-space calibration. No q=0 or diffraction indices assigned. Bootstrap quantifies spatial sampling, not independently repeated acquisitions.\n');
bisb_slot_summary(out);
bisb_slot_anchor_analysis(out);
if isfield(cfg,'geometry'), bisb_slot_geometry(out); end
fprintf('RESULT_DIRECTORY=%s\n',out);
end
