function b1_adaptq_look(run_root, fullq_root)
% B1 adaptive-binning task, stage 1 check: (1) N=3 sums of the aligned native
% channels reproduce the b1_fullq bins; (2) per-channel counts in the B1
% window, to choose the binning; (3) native-resolution stacked spectra and
% q-E map of 590 (no fit, no kinematic correction).
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
out=fullfile(run_root,'stage1'); fig=fullfile(run_root,'figures'); if ~isfolder(fig), mkdir(fig); end
for s={'590_PL2_10w','n0_PL2_10w_repeat'}
 session=s{1}; S=load(fullfile(out,sprintf('channels_%s.mat',session))); C=S.C;
 Y=squeeze(sum(double(S.frames),2)); E=C.E; q=C.q;
 % (1) consistency with b1_fullq
 F=load(fullfile(fullq_root,'stage1',sprintf('frames_%s.mat',session)),'D'); D=F.D;
 for t=[-0.0025 0.0025 0.0100]
  b=find(abs(D.q_center-t)<1e-6); m=ismember(round(q/1e-6),round(D.bins(b).q_members/1e-6));
  assert(isequal(D.E,E),'E axis differs'); d=max(abs(sum(Y(:,m),2)-D.Ysum(:,b)))/max(D.Ysum(:,b));
  fprintf('%s: bin %+.4f vs b1_fullq max rel diff %.1e (%d channels)\n',session,t,d,nnz(m));
 end
 % (2) counts per channel
 w=E>=300&E<=1800; hi=E>=1200&E<=1800; z=abs(E)<=50;
 fprintf('%s per-channel counts (sum over frames):\n  |q|      win300-1800   win1200-1800   zlp(+-50)   counts/4meV@1500(+q)\n',session);
 for qa=[0.0005 0.001 0.0015 0.002 0.003 0.004 0.005 0.0065 0.008 0.01 0.0125 0.015 0.02 0.025 0.03]
  k=find(abs(q-qa)<1e-6); if isempty(k), continue; end
  e15=abs(E-1500)<=2;
  fprintf('  %.4f  %12.3g  %12.3g  %12.3g  %10.3g\n',qa,sum(Y(w,k)),sum(Y(hi,k)),sum(Y(z,k)),mean(Y(e15,k)));
 end
 if ~strcmp(session,'590_PL2_10w'), continue; end
 % (3) stacked native spectra, |q| <= 0.006, both signs
 f=figure('Visible','off','Position',[0 0 1300 650]); tl=tiledlayout(f,1,2);
 sel=find(abs(q)<=0.0060+1e-9); [~,o]=sort(q(sel)); sel=sel(o);
 for side=[-1 1]
  ax=nexttile(tl); hold(ax,'on'); ks=sel(sign(q(sel))==side); [~,o]=sort(abs(q(ks))); ks=ks(o);
  cm=parula(numel(ks)+2); r=E>=250&E<=2000;
  for j=1:numel(ks)
   y=movmean(Y(:,ks(j)),11); y=y/max(y(E>=600&E<=1900)); plot(ax,E(r),y(r)+0.35*(j-1),'Color',cm(j,:),'LineWidth',1);
   text(ax,1960,y(find(E>=1900,1))+0.35*(j-1),sprintf('%.4f',abs(q(ks(j)))),'FontSize',7);
  end
  xlabel(ax,'Energy loss (meV)'); ylabel(ax,'counts (11-pt smoothed, norm. 600-1900) + offset');
  title(ax,sprintf('590, q %s 0, native channels (N=1)',char('<'*(side<0)+'>'*(side>0)))); xlim(ax,[250 2150]); grid(ax,'on');
 end
 exportgraphics(f,fullfile(fig,'native_stack_590_lowq.png'),'Resolution',140); close(f);
 % q-E map, native channels, each column normalised in 400-1900 meV
 f=figure('Visible','off','Position',[0 0 1100 600]); [qs,o]=sort(q); M=movmean(Y(:,o),9,1);
 r=E>=250&E<=2200; M=M(r,:); M=M./max(M(E(r)>=400&E(r)<=1900,:),[],1);
 imagesc(qs,E(r),M); axis xy; colormap(turbo); caxis([0.4 1]); colorbar;
 xlabel('q (A^{-1})'); ylabel('Energy loss (meV)'); title('590 native channels, per-column max-normalised (400-1900 meV), raw counts');
 exportgraphics(f,fullfile(fig,'native_qE_590.png'),'Resolution',140); close(f);
end
end
