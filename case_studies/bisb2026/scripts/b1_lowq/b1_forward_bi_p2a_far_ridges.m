function b1_forward_bi_p2a_far_ridges()
%B1_FORWARD_BI_P2A_FAR_RIDGES  Trace curved detector features outside the P2b window.
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
ids = {'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'};
settings = struct('min_abs_ch',40,'apex_window',[-90 430], ...
    'target',[100 350],'half_width',12,'side_offset',36,'agree_meV',12, ...
    'min_z',3,'slope_bands',[100 160;160 240;240 350]);
failed = {};
for task = ids
    id = task{1};
    t = tic;
    try
        cfg = b1_forward_bi_p2a_config(id);
        local_trace(cfg,settings);
        fprintf('DONE %s far ridges %.2f s\n',id,toc(t));
    catch err
        failed{end+1} = id; %#ok<AGROW>
        fprintf(2,'FAIL %s far ridges\n%s\n',id,getReport(err,'extended'));
    end
end
assert(isempty(failed),'Far-ridge checks failed for %s.',strjoin(failed,', '));
end

function local_trace(cfg,s)
p = cfg.paths;
d = bisb_datasets(p.session);
out = fullfile(cfg.input.output_dir,'b1_forward_bi_far_ridges');
assert(~isfolder(out),'Far-ridge output already exists.');
mkdir(out);
S = load(p.stage1_file,'C','cfg');
C = S.C;
A = qe_align_sequence_frames(d.path,d.zero,d.dq,q_select=Inf, ...
    reference_q=S.cfg.reference_abs_q,alignment_window=S.cfg.zlp_window, ...
    energy_zero="detector_peak",raw_profile_window=S.cfg.raw_profile_window);
assert(isequal(A.alignment.measured_offset_pixels(:),C.offsets(:)));
E = A.E(:); q = A.Q;
k = round(q/d.dq);
odd = mod(C.block,2)==1;
M = squeeze(sum(A.F(:,odd,:),2));
N = squeeze(sum(A.F(:,~odd,:),2));
clear A
J = M+N;
[has,mirror] = ismember(-k,k);
columns = find(abs(k)>=s.min_abs_ch);
rows = {};
for j = columns
    [e,inf] = qe_peak_apex(E,J(:,j),s.apex_window,smooth=5,min_prominence=.1,half_width=12);
    if ~isfinite(e), continue; end
    e1 = qe_peak_apex(E,M(:,j),s.apex_window,smooth=5,min_prominence=.1,half_width=12);
    e2 = qe_peak_apex(E,N(:,j),s.apex_window,smooth=5,min_prominence=.1,half_width=12);
    [c1,b1,z1] = local_contrast(E,M(:,j),e,s);
    [c2,b2,z2] = local_contrast(E,N(:,j),e,s);
    cm = NaN;
    if has(j), cm = local_contrast(E,J(:,mirror(j)),e,s); end
    confirmed = e>=s.target(1) && e<=s.target(2) && abs(e1-e2)<=s.agree_meV && ...
        abs(e-e1)<=s.agree_meV && abs(e-e2)<=s.agree_meV && z1>=s.min_z && z2>=s.min_z;
    rows(end+1,:) = {k(j),q(j),e,e1,e2,c1+c2,b1+b2,(c1+c2)/max(b1+b2,1), ...
        cm,z1,z2,inf.prominence,confirmed}; %#ok<AGROW>
end
T = cell2table(rows,'VariableNames',{'native_k','q','apex_meV','odd_apex_meV','even_apex_meV', ...
    'excess_counts','sideband_counts','excess_over_sideband','mirror_excess_counts','odd_z','even_z','prominence','confirmed'});
writetable(T,fullfile(out,'b1_forward_bi_far_ridge_columns.csv'));
keep = T.confirmed;
for sg = [-1 1]
    m = keep & sign(T.native_k)==sg;
    fprintf('%s side %+d confirmed %d columns (independent half samples)\n',p.session,sg,nnz(m));
    if ~any(m), continue; end
    row = T(m,:);
    gap = [1; find(diff(row.native_k)~=1)+1; height(row)+1];
    [~,seg] = max(diff(gap));
    row = row(gap(seg):gap(seg+1)-1,:);
    fprintf('%s side %+d longest contiguous %d k [%g,%g] q [%.5f,%.5f] E [%.2f,%.2f] counts median %.2f relative %.3f mirror counts %.2f\n', ...
        p.session,sg,height(row),min(row.native_k),max(row.native_k),min(row.q),max(row.q), ...
        min(row.apex_meV),max(row.apex_meV),median(row.excess_counts), ...
        median(row.excess_over_sideband),median(row.mirror_excess_counts,'omitnan'));
    for v = s.slope_bands.'
        take = row.apex_meV>=v(1) & row.apex_meV<=v(2);
        if nnz(take)<3, continue; end
        poly = polyfit(row.native_k(take),row.apex_meV(take),1);
        fprintf('%s side %+d band %g-%g k [%g,%g] slope %.3f meV/ch excess median %.2f ratio %.2f\n', ...
            p.session,sg,v(1),v(2),min(row.native_k(take)),max(row.native_k(take)),poly(1), ...
            median(row.excess_counts(take)),median(row.excess_over_sideband(take)));
    end
end
f = figure('Visible','off','Theme','light','Position',[0 0 1250 600]);
ax = axes(f); m = abs(E)<=360;
imagesc(ax,k,E(m),log10(max(J(m,:),1))); axis(ax,'xy'); hold(ax,'on'); colorbar(ax);
plot(ax,T.native_k(keep),T.apex_meV(keep),'w.','MarkerSize',7);
xlabel(ax,'Signed native q channel'); ylabel(ax,'Energy (meV)');
title(ax,sprintf('%s: off-axis loss maxima confirmed in independent block sums',p.session),'Interpreter','none');
exportgraphics(f,fullfile(out,'b1_forward_bi_far_ridges.png'),'Resolution',120);
close(f);
record = struct('settings',s,'sha256',b1_sha256(s),'code_sha256', ...
    b1_sha256({'file',[mfilename('fullpath') '.m']}),'raw_sha256',C.raw_sha256);
fid = fopen(fullfile(out,'b1_forward_bi_far_ridge_config.json'),'w');
assert(fid>=0); fwrite(fid,jsonencode(record,'PrettyPrint',true)); fclose(fid);
end

function [excess,bg,z] = local_contrast(E,y,e,s)
x = (-s.half_width:4:s.half_width).'+e;
p = sum(interp1(E,y,[x x-s.side_offset x+s.side_offset],'linear'),1);
bg = (p(2)+p(3))/2;
excess = p(1)-bg;
z = excess/sqrt(max(p(1)+(p(2)+p(3))/4,1));
end
