function B = b1_forward_bi_p2c_prepare(P,stage,h)
%B1_FORWARD_BI_P2C_PREPARE  Share the P2b detector operator across matched fits.
p = P.c.paths;
main = load(fullfile(P.archive,sprintf('fit_%s_%s',P.id,P.main),'p2b_results.mat'));
seed = load(fullfile(P.archive,sprintf('fit_%s_%s',P.id,P.seed),'p2b_results.mat'),'res','cfg');
c = seed.cfg;
o = P.o;
c.exit_q = main.cfg.exit_q;
c.variants = c.variants(1);
c.variants.window = [150 1800];
c.variants.empirical_centre = P.empirical;
if strcmp(stage,'odd')
    c.variants.window = [250 1800];
else
    o.n_sy = max(9,1+2*round((P.nsy-1)/2*h/P.o.h));
end
o.h=h;
if P.fastNoise
    N = load(fullfile(p.diag_dir,'noise_model.mat'),'E','g','q_noise');
    [~,ix] = min(abs(c.exit_q(:)-N.q_noise(:).'),[],2);
    c.noise_gain = struct('E',N.E(:),'g',N.g(:,ix));
end
if P.gapless, c.variants.ext = 'gapless'; end
B = b1_p2b_fit([],c,o,struct('prepare',true));
B.seed = seed.res(1).fit;
B.archive = main.mainMu;
B.cfg = c;
B.ocfg = o;
B.stage = stage;
L = load(fullfile(p.input_dir,'p2a_cube.mat'),'cube');
R = load(fullfile(p.response_dir,'response.mat'),'R');
B.cube = L.cube;
B.response = R.R.b6;
fprintf('%s %s operator %.2f s, h=%.6g n_sy=%d\n',P.id,stage,B.build_s,o.h,o.n_sy);
end
