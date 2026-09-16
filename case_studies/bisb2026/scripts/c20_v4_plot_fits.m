function c20_v4_plot_fits(d,path)
f=figure('Visible','off','Position',[100 100 1150 850]); tiledlayout(2,2);
for n=1:2
 z=d.fits(n); nexttile(n); h=plot(z.energy_meV,z.observed,'Color',[.6 .6 .6]); hold on;
 labels="obs";
 if z.success
  h(end+1)=plot(z.energy_meV,z.background,'--','LineWidth',1.1); labels(end+1)="bg";
  for j=1:n, h(end+1)=plot(z.energy_meV,z.components(:,j),'LineWidth',1.2); labels(end+1)="P"+j; end
  h(end+1)=plot(z.energy_meV,z.prediction,'k','LineWidth',1.4); labels(end+1)="total";
 end
 legend(h,labels,'Location','best'); assert(numel(h)==numel(labels)); ylabel('Corrected counts (q mean)');
 title(sprintf('n=%d: %s',n,z.numerical_status),'Interpreter','none');
 nexttile(n+2); plot(z.energy_meV,z.residual); yline(0,'k:'); xlabel('Energy (meV)'); ylabel('Observed - prediction');
end
sgtitle(sprintf('%s | q=%+.4f | native members %s',d.key,d.unit.q_Ainv,mat2str(d.unit.source_channel)),'Interpreter','none');
exportgraphics(f,path); close(f);
end
