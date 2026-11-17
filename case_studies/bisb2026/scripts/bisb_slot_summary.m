function bisb_slot_summary(out)
%BISB_SLOT_SUMMARY Compare strong elastic spots and spatial variation.
a=load(fullfile(out,'analysis.mat')); E=a.E; rng(a.cfg.seed);
ranges=[65 80;89 104;225 232;233 240;369 376;121 200;265 344];
labels={'left_outer';'left_spot';'central_left';'central_right';'right_spot';'between_left';'between_right'};
K=numel(labels); apex=nan(K,1); ci=nan(K,2); success=nan(K,1); bgpeak=nan(K,3); sm=nan(K,3);
f=figure('Visible','off','Color','w','Position',[50 50 1150 650]); tiledlayout(1,2);
for j=1:K
 take=a.qb>=ranges(j,1)&a.qb<=ranges(j,2); X=squeeze(sum(a.S(:,take,:),2)); y=sum(X,2);
 apex(j)=qe_peak_apex(E,y,a.cfg.peak_window,smooth=15,half_width=60);
 for k=1:3
  [bg,r]=qe_powerlaw_tail(E,y,a.cfg.background_windows(k,:));
  if r>0, bgpeak(j,k)=qe_peak_apex(E,y-bg,a.cfg.peak_window,smooth=15,half_width=60); end
  smooths=[11 21 31]; sm(j,k)=qe_peak_apex(E,y,a.cfg.peak_window,smooth=smooths(k),half_width=60);
 end
 rows=squeeze(sum(reshape(X,numel(E),19,19),2)); b=nan(1,a.cfg.bootstrap);
 for k=1:a.cfg.bootstrap
  yy=sum(rows(:,randi(19,1,19)),2); b(k)=qe_peak_apex(E,yy,a.cfg.peak_window,smooth=15,half_width=60);
 end
 ci(j,:)=prctile(b,[2.5 97.5]); success(j)=mean(isfinite(b));
 nexttile(1+(j>5)); hold on; norm=max(movmean(y(E>=400&E<=1200),15)); plot(E/1000,movmean(y,7)/norm,'DisplayName',strrep(labels{j},'_',' '),'LineWidth',1.3);
end
for j=1:2
 nexttile(j); xlim([.15 1.8]); ylim([.3 1.5]); grid on; xlabel('Loss (eV)'); ylabel('Peak-normalized measured intensity'); legend('Location','best');
end
exportgraphics(f,fullfile(out,'06_spot_comparison.png'),'Resolution',180); close(f);
T=table(labels,ranges(:,1),ranges(:,2),apex,ci(:,1),ci(:,2),success,bgpeak(:,1),bgpeak(:,2),bgpeak(:,3),min(sm,[],2),max(sm,[],2), ...
 'VariableNames',{'region','channel_first','channel_last','raw_apex_meV','ci_low','ci_high','success','bg1_meV','bg2_meV','bg3_meV','smooth_min','smooth_max'});
writetable(T,fullfile(out,'spot_comparison.csv'));
take=a.qb>=225&a.qb<=240; X=squeeze(sum(a.S(:,take,:),2)); peak=nan(1,361);
for t=1:361, peak(t)=qe_peak_apex(E,X(:,t),a.cfg.peak_window,smooth=15,half_width=60); end
rat=sum(X(E>=400&E<=800,:),1)./sum(X(E>=1200&E<=1600,:),1);
f=figure('Visible','off','Color','w','Position',[50 50 1000 450]); tiledlayout(1,2);
nexttile; pm=reshape(peak,19,19)'; im=imagesc(pm); im.AlphaData=isfinite(pm); axis image; colorbar; title('Central spot peak apex (meV)'); xlabel('Scan x pixel'); ylabel('Scan y pixel');
nexttile; imagesc(reshape(rat,19,19)'); axis image; colorbar; title('I(0.4-0.8 eV) / I(1.2-1.6 eV)'); xlabel('Scan x pixel'); ylabel('Scan y pixel');
exportgraphics(f,fullfile(out,'07_spatial_check.png'),'Resolution',180); close(f);
fid=fopen(fullfile(out,'spatial_summary.txt'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'Central pixels with detected peak %d/361\nPeak 5/50/95 percentiles (meV): %g %g %g\nBand ratio 5/50/95 percentiles: %g %g %g\n',sum(isfinite(peak)),prctile(peak,[5 50 95]),prctile(rat,[5 50 95]));
[tf,pr]=islocalmax(a.elastic,'MinProminence',max(a.elastic)*.0005); ii=find(tf); [~,order]=sort(pr(ii),'descend');
fprintf(fid,'Strongest elastic local maxima (channels): '); fprintf(fid,'%g ',a.q(ii(order(1:min(6,end))))); fprintf(fid,'\n');
% Paired spatial-row bootstrap keeps left/right acquisition fluctuations paired.
L=squeeze(a.S(:,find(a.qb==228.5,1),:)); H=squeeze(a.S(:,find(a.qb==236.5,1),:));
L=squeeze(sum(reshape(L,numel(E),19,19),2)); H=squeeze(sum(reshape(H,numel(E),19,19),2));
delta=nan(1,a.cfg.bootstrap);
for k=1:a.cfg.bootstrap
 ix=randi(19,1,19);
 delta(k)=qe_peak_apex(E,sum(H(:,ix),2),a.cfg.peak_window,smooth=15,half_width=60)-qe_peak_apex(E,sum(L(:,ix),2),a.cfg.peak_window,smooth=15,half_width=60);
end
fprintf(fid,'Central right minus left raw apex = %.3f meV; paired 95%% CI %.3f %.3f meV\n',apex(4)-apex(3),prctile(delta,[2.5 97.5]));
ratio_rows=mean(reshape(rat,19,19),1);
fprintf(fid,'Central band ratio first/last scan row: %.6f %.6f; relative change %.3f%%\n',ratio_rows(1),ratio_rows(end),100*(ratio_rows(end)/ratio_rows(1)-1));
disp(T);
end
