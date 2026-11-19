function p = b1_forward_20w_inputs(cfg, out)
%B1_FORWARD_20W_INPUTS  Audit the raw detector and reproduce the input diagnostics.
t = tic;
audit = fullfile(out,'b1_forward_20w_detector.csv');
assert(~isfile(audit),'Detector audit already exists.');
p = cfg.paths;
d = bisb_datasets(p.session);
assert(strcmp(b1_sha256({'file',d.path}),d.sha256),'Raw SHA-256 mismatch.');
S = load(p.stage1_file,'C');
C = S.C;
raw = read_npy(d.path);
src = cfg.bad_column + (-3:3);
tol = double(intmax('uint32'));
rows = zeros(numel(src), 12);
for k = 1:numel(src)
    x = double(squeeze(raw(:,src(k),:)));
    valid = isfinite(x) & x < tol;
    counts = sum(x,2,'omitnan');
    countsClean = x;
    countsClean(~valid) = 0;
    byFrame = sum(countsClean,2);
    rows(k,:) = [src(k),nnz(~isfinite(x)),nnz(x>=tol),nnz(x==0),nnz(x<0), ...
        nnz(all(~valid,2)),nnz(sum(valid,2)<size(x,2)), ...
        sum(counts,'omitnan'),sum(byFrame),median(byFrame),min(byFrame),max(byFrame)];
end
T = array2table(rows,'VariableNames',{'column','nonfinite','sentinel','zero','negative', ...
    'all_invalid_frames','affected_frames','raw_sum','valid_sum','median_frame_sum','min_frame_sum','max_frame_sum'});
T.adjacent_valid_sum_ratio = T.valid_sum / mean(T.valid_sum(ismember(T.column,cfg.bad_column+[-1 1])));
writetable(T,audit);
disp(T)
old = load(fullfile(p.diag_dir,'common_ref.mat'),'ref','validq','ho');
validq = true(1,size(raw,2));
for k = 1:numel(validq)
    x = raw(:,k,:);
    validq(k) = all(isfinite(x(:)) & x(:)<tol);
end
fprintf('Raw invalid channels = %s; archive invalid channels = %s\n', ...
    mat2str(find(~validq)),mat2str(setdiff(1:numel(validq),old.validq)));
assert(~any(ismember(C.source_channel,find(~validq))),'Invalid channel in archive stage1.');
L = load(p.l1_file);
validq = validq & all(isfinite(L.d.qe.intensity),1);
assert(isequal(find(validq),old.validq(:).'),'Archive validity differs from raw and L1.');
F = load(p.stage1_file,'frames');
ref = zeros(C.T,1);
hidx = C.support(C.E>=old.ho(1)&C.E<=old.ho(2));
for k = 1:C.T
    f = double(squeeze(raw(k,validq,:)));
    ref(k) = sum(f(:))-sum(sum(f(:,hidx+C.offsets(k))));
    x = reshape(squeeze(raw(k,C.source_channel,C.support+C.offsets(k))).',size(F.frames,1),1,[]);
    assert(isequal(x,F.frames(:,k,:)),'Archive stage1 frame %d differs from raw.',k);
end
refError = max(abs(ref-old.ref(:)));
assert(refError<=8*eps(max(ref)),'Archive common reference does not match raw.');
clear raw L
fprintf('Verified %d aligned frames; reference max diff %.3g; raw SHA-256 %s\n',C.T,refError,d.sha256);
i = b1_p2a_inputs_config(p);
i.version = [i.version '+' cfg.version];
p.input_dir = b1_p2a_inputs([],i);
q = reshape([-i.key_q;i.key_q],1,[]);
E = C.E;
g = zeros(numel(E),numel(q));
sigma = g;
for k = 1:numel(q)
    channel = abs(C.q-q(k))<1e-10;
    [sigma(:,k),g(:,k)] = b1_p1b_noise(F.frames(:,:,channel),E,ref,old.ho,'first_diff',21);
end
save(fullfile(p.input_dir,'b1_forward_20w_noise.mat'),'E','q','g','sigma','ref');
clear F
r = b1_p2a_response_config(p);
p.response_dir = b1_p2a_response([],r);
o = b1_p2a_operator_config(p);
p.operator_dir = b1_p2a_operator([],[],o);
base = b1_forward_config();
M = [b1_forward_20w_input_metrics(base);b1_forward_20w_input_metrics(p)];
writetable(M,fullfile(p.input_dir,'b1_forward_20w_vs_590_inputs.csv'));
disp(M);
fprintf('S1 input %s\nS1 response %s\nS1 operator %s\nS1 completed %.2f s\n', ...
    p.input_dir,p.response_dir,p.operator_dir,toc(t));
save(fullfile(out,'b1_forward_20w_inputs.mat'),'p','i','r','o');
end
