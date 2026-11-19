function b1_forward_bi_p2c_low(P,out,phase,deadline)
%B1_FORWARD_BI_P2C_LOW  Run low-window fits, background controls, and h scan.
t0=tic;
fid=fopen(fullfile(P.run,'step1_choice.csv'),'r');
assert(fid>0,'Missing Step 1 model choice.');
Dec=textscan(fid,'%q%q%q','Delimiter',',','HeaderLines',1);
fclose(fid);
row=find(strcmp(Dec{1},P.id));
assert(isscalar(row),'Expected one Step 1 model choice for %s.',P.id);
P.maxIter=60;
P.gain=Dec{2}{row};
P.odd_source=Dec{3}{row};
P.qfloor=0;
if strcmp(P.id,'Bi1117_10w')
    % P2a gain-proxy stable wings: avoid the elastic-centre-dominated pair.
    P.qfloor=2*bisb_datasets(P.id).dq;
end
kind=regexp(phase,'^low_(pilot|primary|controls|hscan|threshold|z_nested)','match','once');
if endsWith(phase,'_followup'), P.skip_checked_nested=true; end
V=local_tasks(P,kind);
nFail=0;
B=[];
lastKey='';
for k=1:size(V,1)
    if datetime('now','TimeZone','UTC')>=deadline
        fprintf(2,'Wall-clock cutoff %s. Unstarted: %s\n',P.id,strjoin(V(k:end,1),', '));
        break
    end
    try
        t=tic;
        name=V{k,1}; p=V{k,2}; h=V{k,3};
        folder=fullfile(out,name);
        assert(~isfolder(folder),'Results already exist for %s',folder);
        mkdir(folder);
        key={h,p.seed,p.empirical,p.fastNoise,p.gapless};
        if ~isequal(key,lastKey)
            B=[];
            B=b1_forward_bi_p2c_prepare(p,'low',h);
            lastKey=key;
        end
        % P2b and odd-fold fits provide starting points only. The low
        % fit takes the same archived material seed for all families.
        tasks=V{k,4};
        fails=local_folds(B,p,tasks,folder,name);
        nFail=nFail+fails;
        fprintf('DONE %s %s in %.2f s (%d failures)\n',P.id,name,toc(t),fails);
        if startsWith(phase,'low_pilot')
            fprintf('Low-window pilot %.2f s. 6 background variants + 4--5 h points + up to 6 controls per group: %.2f h at pilot rate over 2 workers, plus full-grid and second-fold overhead.\n', ...
                toc(t),4*17*toc(t)/7200);
        end
    catch err
        nFail=nFail+1;
        fprintf(2,'FAILED %s %s: %s\n',P.id,V{k,1},getReport(err,'extended','hyperlinks','off'));
    end
end
fprintf('%s %s %.2f s, %d failures\n',P.id,phase,toc(t0),nFail);
assert(nFail==0,'%d low-window task(s) failed.',nFail);
end

function V=local_tasks(P,phase)
main={6,true;1,false};
V=cell(0,4);
switch phase
    case 'low_pilot'
        P.briefGrid=true;
        V={'main',P,P.o.h,{6,true}};
    case 'low_primary'
        V={'main',P,P.o.h,main};
        for x={'power1','power2','power3','anchor'}
            p=P; p.tail=x{1}; p.briefGrid=true;
            V(end+1,:)={x{1},p,P.o.h,main}; %#ok<AGROW>
        end
        p=P; p.phonon=false; p.briefGrid=true;
        V(end+1,:)={'no_phonon',p,P.o.h,main};
        for x={'power1','anchor'}
            p.tail=x{1};
            V(end+1,:)={['no_phonon_' x{1}],p,P.o.h,{6,true}}; %#ok<AGROW>
        end
    case 'low_controls'
        P.briefGrid=true;
        p=P; p.window=P.control;
        V={'edge',p,P.o.h,main};
        if strcmp(P.id,'Bi1107_20w')
            p=P; p.qmax=max(abs(P.c.exit_q)); p.seed=P.main;
            V(end+1,:)={'full_q',p,P.o.h,main};
            p.window=P.control;
            V(end+1,:)={'full_q_edge',p,P.o.h,main};
            p=P; p.seed='b8a5cced';
            V(end+1,:)={'one_branch',p,P.o.h,main};
        elseif strcmp(P.id,'Bi1117_10w')
            p=P; p.qfloor=0;
            V(end+1,:)={'all_pairs_150',p,P.o.h,main};
        elseif strcmp(P.id,'Bi1117_20w')
            p=P; p.empirical=true; p.seed='ee5f4271';
            V(end+1,:)={'empirical_center',p,P.o.h,main};
        end
        p=P; p.gain='none';
        V(end+1,:)={'symmetric',p,P.o.h,{6,true}};
        if strcmp(P.gain,'gainE')
            p=P; p.gain='shift';
            V(end+1,:)={'shift_proxy',p,P.o.h,{6,true}};
        end
        p=P; p.fastNoise=true;
        V(end+1,:)={'fast_noise',p,P.o.h,{6,true}};
        p=P; p.gapless=true;
        V(end+1,:)={'gapless_extension',p,P.o.h,{6,true}};
        p=P; p.fullB1=true;
        V(end+1,:)={'free_B1',p,P.o.h,{6,true}};
    case 'low_hscan'
        P.briefGrid=true;
        inv=readtable(fullfile(P.archive,'compare_16faaaff','fit_inventory.csv'),'TextType','string');
        for h=P.hScan(abs(P.hScan-P.o.h)>1e-9)
            p=P;
            t=inv.id==string(P.id) & abs(inv.h-h)<1e-9 & startsWith(inv.label,"hscan") & inv.variant=="main";
            assert(nnz(t)==1,'Expected one archived mixture-y P2b fit at h=%g.',h);
            sha=char(inv.sha256(t));
            p.seed=sha(1:8);
            V(end+1,:)={sprintf('h_%.6f',h),p,h,{6,true}}; %#ok<AGROW>
        end
    case 'low_threshold'
        s=readtable(fullfile(P.run,'low_primary_profile',P.id,'main','scores_6.csv'), ...
            'TextType','string','Delimiter',',');
        e=s.parameter(s.family=="interband");
        assert(isscalar(e) && isfinite(e),'Expected one nominal interband threshold.');
        P.thresholdGrid=unique([P.thresholdGrid,max(0,e-50):5:min(800,e+50)]);
        V={'threshold_refined',P,P.o.h,{6,true}};
    case 'low_z_nested'
        V=local_nested_tasks(P);
    otherwise
        error('Unknown low-window batch %s.',phase);
end
end

function V=local_nested_tasks(P)
% Retest every nominal comparison, plus detected local-minimum failures.
V=cell(0,4);
source={'low_primary_profile','low_controls','low_hscan','low_threshold'};
kind={'low_primary','low_controls','low_hscan','low_threshold'};
for k=1:numel(kind)
    T=local_tasks(P,kind{k});
    for j=1:size(T,1)
        name=T{j,1};
        path=fullfile(P.run,source{k},P.id,name);
        if isfield(P,'skip_checked_nested') && P.skip_checked_nested
            checked=fullfile(P.run,'low_z_nested',P.id,name);
            held=T{j,4}(:,1);
            complete=cellfun(@(b) isfile(fullfile(checked,sprintf('scores_%d.csv',b))),held);
            if all(complete), continue; end
        end
        take=strcmp(name,'main');
        for held=[6 1]
            fname=fullfile(path,sprintf('scores_%d.csv',held));
            if ~isfile(fname), continue; end
            S=readtable(fname,'TextType','string','Delimiter',',');
            none=S.objective(S.family=="none");
            if any(S.family~="none" & S.W==0 & S.delta_objective>max(.25,1e-5*none))
                take=true;
            end
        end
        if take
            task=T(j,:);
            task{2}.nesting_check=true;
            V(end+1,:)=task; %#ok<AGROW>
        end
    end
end
end

function nFail=local_folds(B,P,tasks,out,name)
nFail=0;
for j=1:size(tasks,1)
    held=tasks{j,1};
    doLimits=tasks{j,2};
    try
        t=tic;
        seed=B.seed;
        if ~strcmp(P.gain,'none')
            s=load(fullfile(P.odd_source,sprintf('fit_%s_block%d.mat',P.gain,held)),'F');
            P.fixedGain=s.F.gain;
            if strcmp(P.seed,P.main), seed=s.F; end
        end
        fprintf('START %s %s held %d, q %.5g--%.5g, E %d--%d\n',P.id,name,held,P.qfloor,P.qmax,P.window);
        [M,prof,limits]=b1_forward_bi_p2c_search(B,P,held,seed,doLimits);
        score=local_scores(M,B,P);
        writetable(score,fullfile(out,sprintf('scores_%d.csv',held)));
        writetable(prof,fullfile(out,sprintf('profile_%d.csv',held)));
        writetable(limits,fullfile(out,sprintf('limits_%d.csv',held)));
        tag=b1_sha256(struct('fit',P,'geometry',B.ocfg,'held',held));
        meta=struct('config',P,'config_sha256',tag,'operator',B.ocfg,'baseline_head', ...
            strtrim(fileread(fullfile(P.run,'git_head_before.txt'))),'elapsed_s',toc(t), ...
            'geometry_s',B.build_s,'node_grid',B.cfg.nodes,'held',held, ...
            'fit_window_meV',P.window,'ny',B.ocfg.n_sy,'code',local_sha());
        fid=fopen(fullfile(out,sprintf('config_%d.json',held)),'w');
        assert(fid>0); fwrite(fid,jsonencode(meta,'PrettyPrint',true)); fclose(fid);
        if held==6
            Data=struct('E',B.E,'q',B.q,'Y',B.Y,'S2',B.S2,'exposure',B.exposure);
            save(fullfile(out,'predictions_6.mat'),'Data','M','P');
            if strcmp(name,'main'), b1_forward_bi_p2c_figures(P,Data,M,out); end
        end
        fprintf('DONE %s %s block %d %.2f s %s\n',P.id,name,held,toc(t),tag);
    catch err
        nFail=nFail+1;
        fprintf(2,'FAILED %s %s block %d: %s\n',P.id,name,held,getReport(err,'extended','hyperlinks','off'));
    end
end
end

function T=local_scores(M,B,P)
names={'none','interband','drude'};
r={};
for k=1:numel(names)
    F=M{k};
    pars=[NaN NaN]; mf=NaN; amp=NaN;
    if k>1, pars=F.shape; mf=F.commonFrac; amp=F.amplitude; end
    r(k,:)={names{k},F.held,P.window(1),P.window(2),nnz(any(any(F.mask,1),3)),B.ocfg.h, ...
        F.chi2,M{1}.chi2-F.chi2,F.objective,M{1}.objective-F.objective,F.chi2dof,F.dof,F.n_held, ...
        F.held_chi2n,F.train_valley,F.held_valley,F.W,F.relativeW,F.W*mf,F.Wlow,pars(1),pars(2),amp,F.flag,F.iterations, ...
        mat2str(F.tailCoef,5),mat2str(F.phCoef,5)}; %#ok<AGROW>
end
T=cell2table(r,'VariableNames',{'family','held','Elo','Ehi','n_q','h','train_chi2','train_delta_chi2', ...
    'objective','delta_objective','train_chi2dof','dof','n_held','held_chi2n','train_valley_chi2n','held_valley_chi2n', ...
    'W','W_over_low','W_common','Wlow','parameter','width','amplitude','exitflag','iterations','tail_scales','phonon_scales'});
end

function S=local_sha()
names={'config','prepare','gain','background','material','fit','search','low','figures'};
S=struct();
for k=1:numel(names)
    n=names{k};
    S.(n)=b1_sha256({'file',which(['b1_forward_bi_p2c_' n])});
end
S.p2b_fit=b1_sha256({'file',which('b1_p2b_fit')});
end
