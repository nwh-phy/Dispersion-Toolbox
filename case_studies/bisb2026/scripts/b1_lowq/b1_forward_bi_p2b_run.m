function b1_forward_bi_p2b_run(phase,ids)
%B1_FORWARD_BI_P2B_RUN  Run P2b fit tasks serially with independent logs.
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
if nargin<2, ids = {'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'}; end
tasks = cell(0,3);
for k = 1:numel(ids)
    [c,o,hScan] = b1_forward_bi_p2b_config(ids{k});
    tasks = [tasks;local_tasks(phase,c,o,hScan)]; %#ok<AGROW>
end
clock = tic;
failed = 0;
started = strtrim(fileread(fullfile(tasks{1,2}.paths.p2b_dir,'started_utc.txt')));
deadline = datetime(started,'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''','TimeZone','UTC')+hours(4)+minutes(40);
for k = 1:size(tasks,1)
    if datetime('now','TimeZone','UTC')>=deadline
        fprintf(2,'Wall-clock cutoff. Unstarted tasks: %s\n',strjoin(tasks(k:end,1),', '));
        break
    end
    try
        c = tasks{k,2};
        o = tasks{k,3};
        [c,expected] = local_unfinished(c);
        if isfile(fullfile(expected,'config.json'))
            fprintf('Existing complete fit %s: %s\n',tasks{k,1},expected);
            continue
        end
        tk = tic;
        fprintf('START %s %s h=%.6g n_sy=%d\n',c.paths.session,tasks{k,1},o.h,o.n_sy);
        result = b1_p2b_fit([],c,o);
        elapsed = toc(tk);
        fprintf('DONE %s %s %.2f s %s\n',c.paths.session,tasks{k,1},elapsed,result);
        if strcmp(phase,'pilot')
            fprintf('Pilot %.2f s. Tier 1 has 49 outputs and 52 optimizer starts. At pilot rate %.2f h serial before geometry and postprocessing differences.\n', ...
                elapsed,52*elapsed/3600);
        end
    catch err
        failed = failed+1;
        fprintf(2,'FAILED %s: %s\n',tasks{k,1},getReport(err,'extended','hyperlinks','off'));
    end
end
fprintf('%s total %.2f s, %d failures\n',phase,toc(clock),failed);
assert(failed==0,'One or more P2b fit tasks failed.');
end

function [c,folder] = local_unfinished(c)
version = c.version;
for k = 0:99
    if k>0, c.version = sprintf('%s+retry%d',version,k); end
    tag = b1_sha256(c);
    folder = fullfile(c.paths.p2b_dir,[c.paths.fit_prefix tag(1:8)]);
    if ~isfolder(folder) || isfile(fullfile(folder,'config.json')), return; end
end
error('b1_forward_bi_p2b_run:Retries','Too many incomplete directories for %s.',version);
end

function T = local_tasks(phase,c,o,hScan)
starts = c.starts;
main = c.variants(1);
T = cell(0,3);
switch phase
    case 'pilot'
        c.variants = main;
        c.starts = starts(1);
        T = {'main',c,o};
    case 'tier1'
        names = {'main','y: Gaussian','broad component off for B1','kernel: exit','low-Q extension: gapless'};
        for k = 1:numel(names)
            v = c;
            v.variants = c.variants(strcmp({c.variants.name},names{k}));
            v.starts = starts(1);
            if k>1, v.version = [c.version '+' names{k}]; end
            T(end+1,:) = {names{k},v,o}; %#ok<AGROW>
        end
        for k = 2:numel(starts)
            v = c;
            v.variants = main;
            v.starts = starts(k);
            v.version = [c.version '+start-' starts(k).name];
            T(end+1,:) = {starts(k).name,v,o}; %#ok<AGROW>
        end
        if strcmp(c.paths.session,'Bi1117_20w')
            v = c;
            v.variants = main;
            v.variants.name = 'kernel: empirical centre';
            v.variants.empirical_centre = true;
            v.starts = starts(1);
            v.version = [c.version '+empirical'];
            T(end+1,:) = {'empirical',v,o};
        end
        if strcmp(c.paths.session,'Bi1107_20w')
            v = c;
            v.variants = main;
            v.variants.name = 'one branch';
            v.n_branch = 1;
            v.bounds_E0 = [150 2500];
            v.starts = struct('name',{'single low','single high'}, ...
                'E0',{starts(1).E0(1,:),starts(1).E0(2,:)}, ...
                'G',{starts(1).G(1,:),starts(1).G(2,:)});
            v.version = [c.version '+single'];
            T(end+1,:) = {'single',v,o};
            v = c;
            v.variants = main;
            v.starts = starts;
            v.exit_q = c.exit_q(abs(c.exit_q)<=12*c.nodes.E0(1)+1e-10);
            v.version = [c.version '+12channels'];
            T(end+1,:) = {'12channels',v,o};
        end
        for h = hScan(abs(hScan-o.h)>1e-10)
            [v,w] = b1_p2b_slot_config('hscan',h,c,o);
            v.variants = v.variants(1);
            v.starts = starts(1);
            T(end+1,:) = {'hscan',v,w}; %#ok<AGROW>
        end
    case 'tier2'
        g = c.variants(strcmp({c.variants.name},'y: Gaussian'));
        for h = hScan(abs(hScan-o.h)>1e-10)
            [v,w] = b1_p2b_slot_config('hscan',h,c,o);
            v.variants = g;
            v.starts = starts(1);
            v.version = [v.version '+Gaussian'];
            T(end+1,:) = {'hscan Gaussian',v,w}; %#ok<AGROW>
        end
        names = {c.variants.name};
        variants = c.variants(ismember(names,{'window 300-1800, no tail','upper edge 1700','kernel: integrated'}) | startsWith(names,'+/-'));
        pearson = main;
        pearson.ytype = 'narrow_pearson';
        pearson.name = 'y: narrow Pearson';
        for V = [variants pearson]
            v = c;
            v.variants = V;
            v.starts = starts(1);
            v.version = [c.version '+' V.name];
            T(end+1,:) = {V.name,v,o}; %#ok<AGROW>
        end
        v = c;
        v.variants = main;
        v.starts = starts(1);
        v.lobo_block = 1;
        v.version = [c.version '+held-1'];
        T(end+1,:) = {'held-1',v,o};
        v = c;
        v.variants = main;
        v.starts = starts(1);
        last = .006;
        if c.nodes.E0(1)>=.0005, last = .008; end
        for n = {'E0','G','A'}, v.nodes.(n{1})(end+1)=last; end
        for n = {'E0','G'}, v.starts.(n{1})(:,end+1)=v.starts.(n{1})(:,end); end
        v.Q_curve = unique([v.Q_curve v.nodes.E0]);
        v.version = [c.version '+outer-node'];
        T(end+1,:) = {'outer-node',v,o};
        v = c;
        v.variants = main;
        v.starts = starts(1);
        N = load(fullfile(c.paths.diag_dir,'noise_model.mat'),'E','g','q_noise');
        [~,near] = min(abs(c.exit_q(:)-N.q_noise(:).'),[],2);
        v.noise_gain = struct('E',N.E(:),'g',N.g(:,near));
        v.version = [c.version '+fast-noise'];
        T(end+1,:) = {'fast-noise',v,o};
        if strcmp(c.paths.session,'Bi1117_20w')
            v = b1_p2b_slot_config('hscan',.0066,c,o);
            v.variants = main;
            v.starts = starts(1);
            [~,folder] = local_unfinished(v);
            p = fullfile(folder,'p2b_results.mat');
            S = load(p,'res');
            e = S.res(1).fit.E0(1,c.nodes.E0>=.001);
            if max(e)-min(e)<=50
                [v,w] = b1_p2b_slot_config('hscan',.010,c,o);
                v.variants = main;
                v.starts = starts(1);
                T(end+1,:) = {'hscan 0.010',v,w};
            else
                fprintf('h=0.010 conditional scan skipped: E0 span at 0.0066 is %.1f meV (threshold 50).\n',max(e)-min(e));
            end
        end
    otherwise
        error('b1_forward_bi_p2b_run:Phase','Unknown phase %s.',phase);
end
end
