function out = b1_p2a_response(input_dir)
% P2a step 2.1: slice response from the block count cube (p2a_cube.mat of the
% corrected input run). For every slice r (native channel, |q| <= 0.0065) and
% block b the elastic line is fitted with two asymmetric Pearson VII sharing a
% centre (qe_elastic_line_model; window c + [-60, 20] meV, Poisson weights). The
% model gives the slice weight pi_br (area within |E - c| <= 100 meV), the
% unit-area energy kernel R_br and its mode c_br. The core-window counts of step 1
% are divided by the model capture fractions (w = 8, 16) as a check. A joint
% two-Pearson fit to pi_br (common shape, block positions) gives the sub-slice
% density and the tails for the operator. The angle-integrated line is rebuilt
% from the slice models. Six and twelve blocks.
% Output: <input_dir>/response_<config hash>/.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
cfg=b1_p2a_response_config(); tag=b1_sha256(cfg); t0=tic;
out=fullfile(input_dir,sprintf('response_%s',tag(1:8))); if ~isfolder(out), mkdir(out); end
L=load(fullfile(input_dir,'p2a_cube.mat'),'cube'); cube=L.cube; E=cube.E(:); q=cube.q(:).'; dq=cube.dq; nk=numel(q); dE=E(2)-E(1);
S=load(fullfile(input_dir,'inputs_p2a.mat'),'P','wl'); T=size(S.P,1);
opts=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-12,'StepTolerance',1e-10, ...
 'MaxIterations',800,'MaxFunctionEvaluations',30000);
R=struct(); rows={};
for nbk=cfg.blockings
 J=cube.(sprintf('b%d',nbk)).J; bk=ceil((1:T).'/(T/nbk));
 th=cell(nbk,nk); area=zeros(nbk,nk); cen=area; fw=area; tailw=area; chi=area; bnd=false(nbk,nk);
 capt=zeros(nbk,nk,numel(cfg.capture_w)); ratio=capt; area2=area;
 for b=1:nbk
  for r=1:nk
   [t,inf]=qe_elastic_line_model('fit',E,J(:,r,b),cfg.fit_window,opts); th{b,r}=t;
   % model values are counts per energy channel: divide the integral by dE to get counts
   area(b,r)=qe_elastic_line_model('area',t,[-cfg.support cfg.support])/dE;
   area2(b,r)=qe_elastic_line_model('area',t,[-cfg.support_check(2) cfg.support_check(2)])/dE;
   [cen(b,r),fw(b,r)]=local_mode_fwhm(t);
   tailw(b,r)=1-qe_elastic_line_model('area',t,[-cfg.tail_report cfg.tail_report])/dE/area(b,r);
   chi(b,r)=inf.chi2dof; bnd(b,r)=inf.at_bound;
   for j=1:numel(cfg.capture_w)
    w=cfg.capture_w(j); fm=qe_elastic_line_model('area',t,[cen(b,r)-t.c-w cen(b,r)-t.c+w])/dE/area(b,r); capt(b,r,j)=fm;
    core=sum(S.P(bk==b,r,S.wl==w),1); ratio(b,r,j)=core/(fm*area(b,r));
   end
   rows(end+1,:)={nbk,b,q(r),area(b,r),cen(b,r),fw(b,r),tailw(b,r),area2(b,r)/area(b,r)-1,capt(b,r,1),capt(b,r,2), ...
    ratio(b,r,1),ratio(b,r,2),chi(b,r),bnd(b,r)}; %#ok<AGROW>
  end
 end
 R.(sprintf('b%d',nbk))=struct('th',{th},'area',area,'centre',cen,'fwhm',fw,'tail_beyond_60',tailw,'capture',capt, ...
  'core_over_model',ratio,'chi2dof',chi,'at_bound',bnd,'blk',bk,'d',cube.(sprintf('b%d',nbk)).d);
end
ST=cell2table(rows,'VariableNames',{'n_blocks','block','q','pi_counts','centre_meV','fwhm_meV','mass_beyond_40meV', ...
 'extra_mass_60_to_100meV','capture_w8','capture_w16','core_w8_over_model','core_w16_over_model','chi2dof_window','at_bound'});
writetable(ST,fullfile(out,'response_slices.csv'));

% joint two-Pearson fit to the slice weights (six blocks): common shape, block positions
opts2=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-14,'StepTolerance',1e-12,'MaxIterations',800,'MaxFunctionEvaluations',40000);
pin=R.b6.area; nb=size(pin,1); jb=cube.jb(1); u0=jb.shape; for b=1:nb, u0=[u0 jb.logA(b)+log(mean(pin(b,:))/mean(cube.b6.core_profile_w8(b,:))) jb.pos(b)/dq]; end %#ok<AGROW>
lbd=[log(1e-4) -8 log([0.05 0.05 0.55 1 1 0.55]) repmat([-Inf -3],1,nb)]; ubd=[log(1) 8 log([5 5 20 40 40 10]) repmat([Inf 3],1,nb)];
u=lsqnonlin(@(u) local_joint_res(u,pin,q,dq,nb),u0,lbd,ubd,opts2); rj=local_joint_res(u,pin,q,dq,nb);
J6=struct('shape',u(1:8),'logA',u(9:2:end),'pos',u(10:2:end)*dq,'logrms',sqrt(mean(rj.^2,1)),'note', ...
 'two-Pearson fit to the capture-corrected slice weights; p = [logA pos/dq shape] for qe_angular_profile_model');
% same fit to the core-window counts (w = 8) for comparison
pc=cube.b6.core_profile_w8; u0c=jb.shape; for b=1:nb, u0c=[u0c jb.logA(b) jb.pos(b)/dq]; end %#ok<AGROW>
uc=lsqnonlin(@(u) local_joint_res(u,pc,q,dq,nb),u0c,lbd,ubd,opts2);
J6c=struct('shape',uc(1:8),'logA',uc(9:2:end),'pos',uc(10:2:end)*dq,'logrms',sqrt(mean(local_joint_res(uc,pc,q,dq,nb).^2,1)));
jrow={};
for v=1:2
 if v==1, x=J6; nm='capture-corrected slice weights'; else, x=J6c; nm='core-window counts w=8 (step 1)'; end
 st=local_stats([0 0 x.shape],dq);
 for b=1:nb
  jrow(end+1,:)={nm,b,x.pos(b),exp(x.shape(3))*dq,exp(x.shape(4))*dq,2*exp(x.shape(5)),x.shape(2)*dq,exp(x.shape(6))*dq, ...
   exp(x.shape(7))*dq,2*exp(x.shape(8)),st(1),st(2),x.logrms(b)}; %#ok<AGROW>
 end
end
writetable(cell2table(jrow,'VariableNames',{'data','block','narrow_mode_Ainv','narrow_scale_minus','narrow_scale_plus','narrow_2m', ...
 'broad_offset_Ainv','broad_scale_minus','broad_scale_plus','broad_2m','narrow_fwhm_Ainv','broad_fraction','logrms'}),fullfile(out,'response_joint_fit.csv'));

% figures
cb=turbo(nb+2); ex=[find(abs(q)<1e-9) find(abs(q+dq)<1e-9) find(abs(q-dq)<1e-9) find(abs(q-2*dq)<1e-9) find(abs(q+6*dq)<1e-9) find(abs(q-6*dq)<1e-9)];
f=figure('Visible','off','Theme','light','Position',[0 0 1800 1000]); tl=tiledlayout(f,2,3,'TileSpacing','compact');
for i=1:numel(ex)
 r=ex(i); ax=nexttile(tl); hold(ax,'on'); m=E>=-150&E<=150;
 for b=[1 3 6]
  y=cube.b6.J(m,r,b); mdl=qe_elastic_line_model('eval',E(m),R.b6.th{b,r});
  semilogy(ax,E(m),y,'.','Color',cb(b+1,:),'MarkerSize',6,'DisplayName',sprintf('block %d data',b));
  semilogy(ax,E(m),mdl,'-','Color',cb(b+1,:),'DisplayName',sprintf('block %d model',b));
 end
 xline(ax,R.b6.centre(1,r)+cfg.fit_window,'k:'); set(ax,'YScale','log'); grid(ax,'on'); xlim(ax,[-150 150]);
 yl=ylim(ax); ylim(ax,[max(yl(1),max(cube.b6.J(m,r,1))*1e-5) yl(2)*1.5]);
 xlabel(ax,'E (meV)'); ylabel(ax,'counts per 4 meV'); title(ax,sprintf('slice q = %+.4f (fit window dotted)',q(r)));
 if i==1, legend(ax,'FontSize',7,'Location','southwest'); end
end
exportgraphics(f,fullfile(out,'fig_slice_lines.png'),'Resolution',110); close(f);

f=figure('Visible','off','Theme','light','Position',[0 0 1800 1000]); tl=tiledlayout(f,2,2,'TileSpacing','compact');
ax=nexttile(tl); hold(ax,'on');
for b=1:nb
 semilogy(ax,q,pin(b,:)/sum(pin(b,:)),'o','Color',cb(b+1,:),'MarkerSize',4);
 semilogy(ax,q,qe_angular_profile_model('channel',q,dq,[J6.logA(b) J6.pos(b)/dq J6.shape])/sum(pin(b,:)),'-','Color',cb(b+1,:));
end
set(ax,'YScale','log'); grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'\pi_{br} / \Sigma_r \pi_{br}');
title(ax,'(a) slice weights (elastic-line model areas) and joint two-Pearson fit');
ax=nexttile(tl); hold(ax,'on');
for b=1:nb, plot(ax,q,squeeze(R.b6.core_over_model(b,:,1)),'o-','Color',cb(b+1,:),'MarkerSize',3); plot(ax,q,squeeze(R.b6.core_over_model(b,:,2)),'s:','Color',cb(b+1,:),'MarkerSize',3); end
grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'core counts / (capture fraction x model area)'); ylim(ax,[0.9 1.1]);
title(ax,'(b) core windows rebuilt from the slice models (o: w = 8, s: w = 16)');
ax=nexttile(tl); hold(ax,'on');
for b=1:nb, plot(ax,q,R.b6.centre(b,:),'o-','Color',cb(b+1,:),'MarkerSize',3); end
grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'slice line mode c_{br} (meV)'); title(ax,'(c) slice energy centres per block');
ax=nexttile(tl); hold(ax,'on');
for b=1:nb, plot(ax,q,squeeze(R.b6.capture(b,:,1)),'o-','Color',cb(b+1,:),'MarkerSize',3); plot(ax,q,R.b6.fwhm(b,:)/40,'^:','Color',cb(b+1,:),'MarkerSize',3); end
grid(ax,'on'); xlabel(ax,'q (A^{-1})'); ylabel(ax,'capture fraction of \pm8 meV (o); FWHM/40 meV (\Delta)');
title(ax,'(d) capture fraction of the w = 8 window and line width');
exportgraphics(f,fullfile(out,'fig_slice_weights.png'),'Resolution',110); close(f);

f=figure('Visible','off','Theme','light','Position',[0 0 1500 650]); tl=tiledlayout(f,1,2,'TileSpacing','compact');
m=E>=-184&E<=300;
for b=[1 6]
 ax=nexttile(tl); y=sum(cube.b6.J(m,:,b),2); mdl=zeros(nnz(m),1);
 for r=1:nk, mdl=mdl+qe_elastic_line_model('eval',E(m),R.b6.th{b,r}).*(abs(E(m)-R.b6.th{b,r}.c)<=cfg.support); end
 semilogy(ax,E(m),y,'k.',E(m),mdl,'r-',E(m),max(y-mdl,1),'b-'); grid(ax,'on');
 legend(ax,{'data, sum over |q| <= 0.0065','sum of slice elastic models','data - models (inelastic + misfit)'},'Location','northeast','FontSize',8);
 xlabel(ax,'E (meV)'); ylabel(ax,'counts per 4 meV'); title(ax,sprintf('block %d: angle-integrated line rebuilt from the slices',b));
end
exportgraphics(f,fullfile(out,'fig_integrated_line.png'),'Resolution',110); close(f);

save(fullfile(out,'response.mat'),'R','J6','J6c','cfg','q','dq','E');
meta=struct('config',cfg,'config_sha256',tag,'code_sha256',b1_sha256({'file',[mfilename('fullpath') '.m']}), ...
 'line_model_sha256',b1_sha256({'file',which('qe_elastic_line_model')}),'profile_model_sha256',b1_sha256({'file',which('qe_angular_profile_model')}), ...
 'input_dir',input_dir,'created',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')),'runtime_s',toc(t0));
fid=fopen(fullfile(out,'config.json'),'w'); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true)); fclose(fid);
fprintf('slice response written to %s (%.0f s)\n',out,toc(t0));
end


function [c, fw] = local_mode_fwhm(t)
x = (t.c - 60:0.05:t.c + 60).'; y = qe_elastic_line_model('eval', x, t); [ym, i] = max(y); c = x(i);
l = find(y(1:i) < ym / 2, 1, 'last'); r = i - 1 + find(y(i:end) < ym / 2, 1, 'first'); fw = NaN;
if ~isempty(l) && ~isempty(r), fw = x(r) - x(l); end
end

function r = local_joint_res(u, Pj, q, dq, nb)
r = zeros(size(Pj, 2), nb);
for b = 1:nb
    r(:, b) = log(qe_angular_profile_model('channel', q(:), dq, [u(7 + 2 * b) u(8 + 2 * b) u(1:8)])) - log(Pj(b, :).');
end
end

function st = local_stats(p, dq)
qf = (-0.5:dq / 40:0.5).'; c = qe_angular_profile_model('components', qf, dq, p); d = sum(c, 2);
h = max(c(:, 1)) / 2; i1 = find(c(:, 1) >= h, 1, 'first'); i2 = find(c(:, 1) >= h, 1, 'last');
st = [qf(i2) - qf(i1), sum(c(:, 2)) / sum(d)];
end
