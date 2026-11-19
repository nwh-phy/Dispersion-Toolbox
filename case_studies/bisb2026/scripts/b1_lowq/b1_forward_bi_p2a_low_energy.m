function b1_forward_bi_p2a_low_energy(cube,R,d,cfg,out)
%B1_FORWARD_BI_P2A_LOW_ENERGY  Gain-side tail proxies and near-130 meV excess.
E = cube.E(:);
q = cube.q;
tail = {};
phonon = {};
fitRows = {};
dose = cube.b6.d(:);
cum = cumsum(dose)-dose/2;
for j = find(abs(q)<=.005+1e-9)
    excess = zeros(6,3);
    sigma = zeros(6,1);
    for b = 1:6
        J = cube.b6.J(:,j,b);
        c = R.b6.centre(b,j);
        x = E-c;
        t = R.b6.th{b,j};
        energies = (72:4:320).';
        % Unscaled empirical gain reflection and three alternative log-power
        % fits on the measured gain side, plus a loss-side anchor sensitivity.
        mirror = interp1(-x,J,energies,'linear');
        est = [mirror zeros(numel(energies),4)];
        for r = 1:size(cfg.tail_gain,1)
            w = -fliplr(cfg.tail_gain(r,:));
            [~,p,A] = qe_powerlaw_tail(-x,J,w);
            est(:,r+1) = A*energies.^(-p);
        end
        use = x>=68 & x<=84;
        y = interp1(-x,J,x(use),'linear');
        ratio = sum(J(use))/sum(y);
        est(:,5) = ratio*mirror;
        % The 60 meV compact kernel is zero here; the formal continuation
        % of the fitted core is supplied solely as a lack-of-fit control.
        y = interp1(x,J,energies,'linear');
        line = qe_elastic_line_model('eval',energies+c,t);
        area = R.b6.area(b,j);
        for r = 1:numel(cfg.eval_energy)
            v = cfg.eval_energy(r);
            m = abs(energies-v)<=8;
            p = mean(est(m,:),1,'omitnan');
            tail(end+1,:) = {q(j),b,v,mean(y(m)),mean(line(m)),p(1),p(2),p(3),p(4),p(5), ...
                min(p),max(p),(max(p)-min(p))/mean(y(m)), ...
                (max(p)-min(p))/area,mean(y(m))/dose(b),min(p)/dose(b),max(p)/dose(b),ratio}; %#ok<AGROW>
        end

        use = x>=cfg.phonon_window(1)&x<=cfg.phonon_window(2);
        side = (x>=cfg.phonon_sidebands(1,1)&x<=cfg.phonon_sidebands(1,2)) | ...
               (x>=cfg.phonon_sidebands(2,1)&x<=cfg.phonon_sidebands(2,2));
        p = polyfit(log(x(side)),log(max(J(side),1)),1);
        bg = exp(polyval(p,log(x(use))));
        alt = (x>=88&x<=112) | (x>=164&x<=204);
        p = polyfit(log(x(alt)),log(max(J(alt),1)),1);
        bg2 = exp(polyval(p,log(x(use))));
        flux = @(a,z) mean(J(x>=a & x<=z));
        c130 = flux(122,138)/sqrt(flux(92,108)*flux(162,178))-1;
        yy = interp1(-x,J,x(use),'linear');
        yg = interp1(-x,J,x(side),'linear');
        p = polyfit(log(x(side)),log(max(yg,1)),1);
        bgG = exp(polyval(p,log(x(use))));
        rawExcess = sum(J(use)-bg);
        rawAlt = sum(J(use)-bg2);
        gainExcess = sum(yy-bgG);
        g = exp(-.5*((x(use)-130)/(max(R.b6.fwhm(b,j),20)/2.355)).^2);
        g = g/sum(g);
        w = 1./max(J(use),1);
        gauss = sum(w.*g.*(J(use)-bg))/sum(w.*g.^2);
        gaussG = sum(w.*g.*(yy-bgG))/sum(w.*g.^2);
        sigma(b) = sqrt(1/sum(w.*g.^2))/dose(b);
        excess(b,:) = [gauss gaussG rawAlt]./dose(b);
        phonon(end+1,:) = {q(j),b,cum(b)/sum(dose),d.frame_s*sum(cube.b6.blk<=b), ...
            c130,rawExcess,rawAlt,gainExcess,rawExcess/sum(bg),gauss,gaussG, ...
            gauss/dose(b),sigma(b),gaussG/dose(b),rawAlt/dose(b),gauss/area}; %#ok<AGROW>
    end
    x = cum/sum(dose);
    XX = [ones(6,1) x];
    [~,upper] = qr(XX,0);
    slopeFactor = norm(upper.'\[0;1]);
    for k = 1:3
        y = excess(:,k);
        p = XX\y;
        slopeSE = norm(y-XX*p)/2*slopeFactor;
        association = corrcoef(x,y);
        fitRows(end+1,:) = {q(j),k,mean(y),min(y),max(y),p(2),slopeSE, ...
            y(6)-y(1),median(sigma),association(1,2)}; %#ok<AGROW>
    end
end
T = cell2table(tail,'VariableNames',{'q','block','E_meV','measured_loss','formal_core_continuation', ...
    'mirror_gain','gain_power_1','gain_power_2','gain_power_3','loss_anchored_mirror','proxy_min', ...
    'proxy_max','envelope_over_loss','envelope_over_elastic_area','loss_per_dose','proxy_min_per_dose','proxy_max_per_dose','loss_anchor_scale'});
P = cell2table(phonon,'VariableNames',{'q','block','cumulative_dose_fraction','end_exposure_s','legacy_c130', ...
    'excess_counts','alt_excess_counts','gain_excess_counts','excess_fraction_bg','gauss_area_counts','gauss_gain_area_counts', ...
    'gauss_per_dose','gauss_shot_sd_per_dose','gain_gauss_per_dose','alt_excess_per_dose','gauss_over_elastic_area'});
D = cell2table(fitRows,'VariableNames',{'q','model','mean_per_dose','min_per_dose','max_per_dose', ...
    'slope_per_full_dose','slope_se_between_blocks','block6_minus_1','median_shot_sd_loss','r_with_dose'});
writetable(T,fullfile(out,'b1_forward_bi_tail.csv'));
writetable(P,fullfile(out,'b1_forward_bi_phonon_blocks.csv'));
writetable(D,fullfile(out,'b1_forward_bi_phonon_dose.csv'));

f = figure('Visible','off','Theme','light','Position',[0 0 1200 660]);
tl = tiledlayout(f,1,2,'TileSpacing','compact');
ax = nexttile(tl); hold(ax,'on');
m = T.q==0;
colors = turbo(numel(cfg.eval_energy)+2);
for k = 1:numel(cfg.eval_energy)
    e = cfg.eval_energy(k);
    u = m & T.E_meV==e;
    c = colors(k+1,:);
    patch(ax,[T.block(u);flipud(T.block(u))], ...
        [T.proxy_min(u)./dose;flipud(T.proxy_max(u)./dose)],c, ...
        'FaceAlpha',.18,'EdgeColor',c,'HandleVisibility','off');
    plot(ax,T.block(u),T.measured_loss(u)./dose,'.-','Color',c, ...
        'DisplayName',sprintf('%g meV',e));
end
set(ax,'YScale','log'); xlabel(ax,'Block'); ylabel(ax,'Counts per 4 meV bin / common dose');
title(ax,'q=0 measured loss and gain-side proxy envelope');
legend(ax,'Location','best');
ax = nexttile(tl); hold(ax,'on');
u = P.q==0;
errorbar(ax,P.cumulative_dose_fraction(u),P.gauss_per_dose(u),P.gauss_shot_sd_per_dose(u),'o-','DisplayName','Loss-side excess');
plot(ax,P.cumulative_dose_fraction(u),P.gain_gauss_per_dose(u),'s-','DisplayName','Gain-side control');
plot(ax,P.cumulative_dose_fraction(u),P.alt_excess_per_dose(u),'^-','DisplayName','Alternate loss sidebands');
yline(ax,0); legend(ax,'Location','best');
xlabel(ax,'Cumulative common dose fraction'); ylabel(ax,'130 meV area / common dose');
exportgraphics(f,fullfile(out,'b1_forward_bi_low_energy.png'),'Resolution',120);
close(f);
end
