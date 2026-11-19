function b1_forward_bi_p2b_heatmap(ax,q,E,x,limits)
%B1_FORWARD_BI_P2B_HEATMAP  Display native q columns, leaving q=0 blank.
dq = min(abs(q));
fullQ = (round(min(q)/dq):round(max(q)/dq))*dq;
[tf,i] = ismembertol(q,fullQ,1e-10,'DataScale',1);
assert(all(tf),'Output channels must lie on the native signed q grid.');
z = nan(numel(E),numel(fullQ));
z(:,i) = x;
im = imagesc(ax,1e3*fullQ,E,z,limits);
im.AlphaData = isfinite(z);
set(ax,'YDir','normal','Color',[.93 .93 .93]);
xlabel(ax,'q (10^{-3} Å^{-1})'); ylabel(ax,'Energy (meV)');
colorbar(ax);
end
