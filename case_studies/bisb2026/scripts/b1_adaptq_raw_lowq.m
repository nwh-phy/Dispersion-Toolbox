function b1_adaptq_raw_lowq(run_root, sessions, qlist)
% Raw native-channel spectra at small q (frame sum, A1-aligned, no fit, no
% kinematic correction), 0-4000 meV, normalised after the post-ZLP minimum,
% with the B1 apex marked. Prints the apex per channel.
if nargin<2, sessions={'no_PL2_20w_2film','590_PL2_10w'}; end
if nargin<3, qlist=[0.00025 0.0005 0.00075 0.0010 0.0015 0.0020 0.0030]; end
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
figd=fullfile(run_root,'figures'); if ~isfolder(figd), mkdir(figd); end
f=figure('Visible','off','Position',[0 0 1700 650]); tl=tiledlayout(f,1,numel(sessions),'TileSpacing','compact');
for s=1:numel(sessions)
 ses=sessions{s}; cf=fullfile(run_root,'stage1',sprintf('channels_%s.mat',ses)); S=load(cf,'C'); C=S.C; mf=matfile(cf); E=C.E;
 ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); k=0; cm=lines(numel(qlist));
 fprintf('%s (dq %.5f): raw apex after the post-ZLP minimum, per side\n',ses,C.dq);
 for qa=qlist
  i=find(abs(abs(C.q)-qa)<1e-7); if isempty(i), continue; end
  y=zeros(numel(E),1); for j=i(:)', y=y+sum(double(mf.frames(:,:,j)),2); end   % both signs summed
  ys=movmean(y,25); r=E>=120&E<=600; [~,im]=min(ys+1e30*~r); emin=E(im);
  w=E>emin&E<=1800; [~,ia]=max(ys.*w); ea=E(ia);
  ap=nan(1,2); for t=1:numel(i), yt=movmean(sum(double(mf.frames(:,:,i(t))),2),25); [~,ib]=max(yt.*w); ap(t)=E(ib); end
  fprintf('  |q| = %.5f: post-ZLP minimum %4.0f meV, apex %4.0f meV (sides %s)\n',qa,emin,ea,mat2str(round(ap)));
  sc=max(ys(E>emin&E<=4000)); k=k+1;
  plot(ax,E,ys/sc+0.9*(k-1),'Color',cm(k,:),'LineWidth',1.1,'DisplayName',sprintf('|q| = %.5f',qa));
  plot(ax,ea,1+0.9*(k-1),'kv','MarkerFaceColor','k','HandleVisibility','off');
 end
 xlim(ax,[0 4000]); ylim(ax,[0 0.9*k+0.6]); xlabel(ax,'energy loss (meV)'); ylabel(ax,'raw counts, \pm q summed, normalised after post-ZLP minimum (+offset)');
 title(ax,sprintf('%s (dq = %.5f A^{-1}/channel); triangles: apex',strrep(ses,'_','\_'),C.dq)); legend(ax,'Location','northeast','FontSize',8);
end
exportgraphics(f,fullfile(figd,'raw_lowq_20w_vs_590.png'),'Resolution',130); close(f);
end
