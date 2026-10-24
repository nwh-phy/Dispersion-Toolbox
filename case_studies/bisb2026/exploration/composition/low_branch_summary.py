"""Low-component evidence, uncertainty, and displays on actual spectra."""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from low_branch_run import read,local_data,BOOT_NAMES,quantities
from models import Operator,solve_amplitudes
from qbin_fit import local_columns
from kinematic_audit import csv_out,save
from io_utils import digest,dump


def interval(values):
    return np.percentile(values,[2.5,97.5])


def summarize(output,ridge):
    cfg=read(output/'signature.json')['protocol'];tasks=read(output/'tasks.json')
    maps=read(output/'mappings.json');fits={};boots={};rows=[];checks=[]
    for t in tasks:
        dest=output/'fits'/t['key'];f=read(dest/'main.json');fits[t['key']]=f
        path=dest/'bootstrap.npz';b=None
        if path.exists():
            z=np.load(path);assert z['names'].tolist()==BOOT_NAMES
            b=z['values'];boots[t['key']]=b
        a=f['A'];s=f['S'];r=dict(key=t['key'],session=t['session'],q_Ainv=t['q_Ainv'],
            low_E0_eV=a['E0_low_eV'],high_E0_eV=a['E0_high_eV'],
            low_loss_apex_eV=a['low_apex_eV'],low_count_apex_eV=a['low_count_apex_eV'],
            high_count_apex_eV=a['high_count_apex_eV'],low_area_fraction=a['low_B1_area_fraction'],
            low_gamma_eV=a['gamma_low_eV'],high_gamma_eV=a['parameters'][6],
            train_gain_fraction=1-a['cost']/s['cost'],test_gain_fraction=1-a['test_b1_msr']/s['test_b1_msr'],
            single_test_msr=s['test_b1_msr'],double_test_msr=a['test_b1_msr'],boundary=';'.join(a['boundary']))
        for name in BOOT_NAMES[:10]:r[name+'_lo']=np.nan;r[name+'_hi']=np.nan
        flags=[];r['bootstrap_success_fraction']=0.;r['bootstrap_boundary_fraction']=np.nan
        r['bootstrap_apex_inversion_fraction']=np.nan
        if b is None:flags.append('bootstrap_pending')
        else:
            valid=b[:,10]>0;v=b[valid]
            r['bootstrap_success_fraction']=float(valid.mean())
            if len(v):
                for j,name in enumerate(BOOT_NAMES[:10]):r[name+'_lo'],r[name+'_hi']=interval(v[:,j])
                r['E0_separation_lo'],r['E0_separation_hi']=interval(v[:,1]-v[:,0])
                r['bootstrap_boundary_fraction']=float(v[:,11].mean())
                r['bootstrap_apex_inversion_fraction']=float(v[:,12].mean())
                if r['low_E0_eV_hi']-r['low_E0_eV_lo']>cfg['stable_E0_interval_meV']/1000:
                    flags.append('E0_interval_wide')
                if r['low_count_apex_eV_hi']-r['low_count_apex_eV_lo']>cfg['stable_count_apex_interval_meV']/1000:
                    flags.append('count_apex_interval_wide')
                if r['low_area_fraction_lo']<cfg['minimum_low_B1_area_fraction']:flags.append('low_component_weak')
                if r['E0_separation_lo']<cfg['minimum_E0_separation_meV']/1000:flags.append('components_merge')
                if r['bootstrap_apex_inversion_fraction']>.1:flags.append('apex_order_ambiguous')
            if valid.mean()<cfg['minimum_bootstrap_success_fraction']:flags.append('optimizer_unstable')
        if not a['success']:flags.append('main_optimizer_failure')
        if a['boundary']:flags.append('parameter_boundary')
        if a['low_count_apex_on_edge']:flags.append('count_apex_censored')
        if a['count_apex_order_inversion']:flags.append('apex_order_inverted')
        # Near-optimal multistart alternatives diagnose extra minima, rather
        # than using agreement of warm-started bootstrap fits as sole evidence.
        from dispersion_models import endpoints
        near=[endpoints(np.array(c['parameters']),'A')[0][0] for c in a['candidates'] if c['cost']<=a['cost']+2.]
        r['near_best_low_E0_span_eV']=float(np.ptp(near)) if near else np.nan
        if r['near_best_low_E0_span_eV']>.15:flags.append('multiple_low_E0_solutions')
        alternate_peaks=[]
        d=local_data(t);op=Operator(d)
        mask=(op.energy>=cfg['fit_window_meV'][0])&(op.energy<=cfg['fit_window_meV'][1])
        for c in a['candidates']:
            if c['cost']>a['cost']+2.:continue
            col,_=local_columns(op,np.array(c['parameters']),'A')
            coeff,_=solve_amplitudes(col,d['frames'][:,:cfg['train_frames']].astype(float).sum(1),d['sigma_train'],mask)
            alternate_peaks.append(quantities(op,dict(model='A',components=col*coeff[:,None,:],
                E0_low_eV=endpoints(np.array(c['parameters']),'A')[0][0],
                E0_high_eV=endpoints(np.array(c['parameters']),'A')[1][0]))['low_count_apex_eV'])
        r['near_best_count_apex_span_eV']=float(np.ptp(alternate_peaks)) if alternate_peaks else np.nan
        if r['near_best_count_apex_span_eV']>.1:flags.append('multiple_count_apex_solutions')
        r['window_checked']=False;r['window_count_shift_eV']=np.nan;r['window_E0_shift_eV']=np.nan
        r['noK_low_E0_eV']=np.nan;r['noK_low_count_apex_eV']=np.nan
        for name in ['window400','noK']:
            path=dest/f'{name}.json'
            if not path.exists():continue
            c=read(path);ca=c['A']
            d=local_data(t);m=(d['energy']>=400)&(d['energy']<=1800);nt=cfg['train_frames']
            y=d['frames'][:,nt:].astype(float).sum(1)
            def msr(x):
                return float(np.mean(((y-x['test_scale']*np.asarray(x['prediction']))/d['sigma_test'])[:,m]**2))
            checks.append(dict(key=t['key'],session=t['session'],q_Ainv=t['q_Ainv'],variant=name,
                low_E0_eV=ca['E0_low_eV'],low_count_apex_eV=ca['low_count_apex_eV'],
                low_fraction=ca['low_B1_area_fraction'],boundary=';'.join(ca['boundary']),
                common400_test_gain=1-msr(ca)/msr(c['S']),
                baseline_common400_test_gain=1-msr(a)/msr(s)))
            if name=='window400':
                r['window_checked']=True
                r['window_count_shift_eV']=ca['low_count_apex_eV']-a['low_count_apex_eV']
                r['window_E0_shift_eV']=ca['E0_low_eV']-a['E0_low_eV']
                if abs(r['window_count_shift_eV'])>.1 or abs(r['window_E0_shift_eV'])>.25:
                    flags.append('fit_window_sensitive')
            else:
                r['noK_low_E0_eV']=ca['E0_low_eV'];r['noK_low_count_apex_eV']=ca['low_count_apex_eV']
        r['conditional_stable']=not flags
        r['predictive_supported']=bool(r['test_gain_fraction_lo']>0)
        r['supported']=bool(r['conditional_stable'] and r['predictive_supported'])
        r['status']='supported' if r['supported'] else ('candidate' if r['conditional_stable'] else 'unresolved')
        # A spectral apex can remain constrained while E0 trades off with
        # damping. Preserve separate quality labels for these observables.
        count_flags=[x for x in flags if x not in ['E0_interval_wide','multiple_low_E0_solutions','fit_window_sensitive']]
        if abs(r['window_count_shift_eV'])>.1:count_flags.append('fit_window_sensitive')
        r['count_conditional_stable']=not count_flags
        r['count_supported']=bool(not count_flags and r['predictive_supported'])
        r['count_status']='supported' if r['count_supported'] else ('candidate' if not count_flags else 'unresolved')
        r['count_flags']=';'.join(count_flags)
        r['flags']=';'.join(flags);rows.append(r)
    # CSV columns must also exist in preliminary summaries with missing bootstraps.
    for r in rows:
        r.setdefault('E0_separation_lo',np.nan);r.setdefault('E0_separation_hi',np.nan)
    mapped=[];lookup={r['key']:r for r in rows}
    for m in maps:mapped.append({**m,**lookup[m['key']]})
    csv_out(output/'low_component_points.csv',mapped)
    csv_out(output/'unique_low_component_points.csv',rows)
    if checks:csv_out(output/'window_and_noK_checks.csv',checks)
    deltas=[]
    for sid in ['590','n0','20w']:
        for sign in [-1,1]:
            for n in [3,5]:
                group=[r for r in mapped if r['session']==sid and r['high_channels']==n and r['q_Ainv']*sign>0]
                a=min(group,key=lambda r:abs(abs(r['q_Ainv'])-.001))
                for target in [.005,.010,.015]:
                    z=min(group,key=lambda r:abs(abs(r['q_Ainv'])-target))
                    for j,name in [(0,'E0'),(3,'count_apex')]:
                        field='low_E0_eV' if j==0 else 'low_count_apex_eV'
                        dr=dict(session=sid,sign=sign,high_channels=n,energy_kind=name,
                            q_start_Ainv=a['q_Ainv'],q_end_Ainv=z['q_Ainv'],
                            start_eV=a[field],end_eV=z[field],delta_eV=z[field]-a[field],
                            start_status=a['status'] if j==0 else a['count_status'],
                            end_status=z['status'] if j==0 else z['count_status'],delta_lo_eV=np.nan,
                            delta_hi_eV=np.nan,bootstrap_pairs=0)
                        if a['key'] in boots and z['key'] in boots:
                            ba=boots[a['key']];bz=boots[z['key']];valid=(ba[:,10]>0)&(bz[:,10]>0)
                            dr['bootstrap_pairs']=int(valid.sum())
                            if valid.any():dr['delta_lo_eV'],dr['delta_hi_eV']=interval(bz[valid,j]-ba[valid,j])
                        deltas.append(dr)
    csv_out(output/'low_component_rise.csv',deltas)
    dump(output/'summary.json',dict(unique_q_fits=len(rows),bootstrap_complete=len(boots),
        supported=sum(r['supported'] for r in rows),candidate=sum(r['status']=='candidate' for r in rows),
        count_supported=sum(r['count_supported'] for r in rows),
        count_candidate=sum(r['count_status']=='candidate' for r in rows),
        unresolved=sum(r['status']=='unresolved' for r in rows),summary_code_sha256=digest(Path(__file__))))
    figures(output,ridge,mapped,fits,tasks,checks)


def figures(output,ridge,rows,fits,tasks,checks):
    font=Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    ridge_rows=list(csv.DictReader((ridge/'peaks.csv').open()))
    def rr(sid,n,sign):
        return sorted([r for r in rows if r['session']==sid and r['high_channels']==n and r['q_Ainv']*sign>0],
                      key=lambda r:abs(r['q_Ainv']))
    for field,title,name,limits in [('low_count_apex_eV','低能成分的计数峰顶',
            '01_low_component_count_dispersion',(.25,1.65)),('low_E0_eV','低能振子的 E₀',
            '02_low_component_E0',(.3,1.85))]:
        fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
        for col,sid in enumerate(['590','n0','20w']):
            for row,sign in enumerate([1,-1]):
                ax=axes[row,col]
                if field=='low_count_apex_eV':
                    total=[r for r in ridge_rows if r['session']==sid and r['high_channels']=='5' and float(r['q_Ainv'])*sign>0]
                    ax.plot([abs(float(r['q_Ainv'])) for r in total],[float(r['peak_meV'])/1000 for r in total],
                            '.',color='#c5c5c5',ms=4,zorder=0)
                for n,color in [(3,'#be793b'),(5,'#187896')]:
                    for r in rr(sid,n,sign):
                        q=abs(r['q_Ainv']);native=q<=.006+1e-10
                        if n==5 and native:continue
                        c='#333333' if native else color;p=r[field];lo=r[field+'_lo'];hi=r[field+'_hi']
                        ax.plot([q,q],[lo,hi],color=c,lw=.7,alpha=.6)
                        status=r['count_status'] if field=='low_count_apex_eV' else r['status']
                        supported=r['count_supported'] if field=='low_count_apex_eV' else r['supported']
                        marker='x' if status=='unresolved' else 'o'
                        ax.plot(q,p,marker,ms=4,color=c,mfc=c if supported else 'white',mew=.9)
                ax.axvspan(0,.006,color='#888888',alpha=.07)
                ax.set(title=f'{sid} | {"正 q" if sign>0 else "负 q"}',xlabel='|q| / Å⁻¹',
                       ylabel=title+' / eV',ylim=limits,xlim=(0,.020))
                ax.grid(alpha=.13)
        handles=[Line2D([],[],marker='o',ls='none',color=c,label=l) for c,l in
                 [('#333333','小 q 单通道'),('#be793b','高 q 合并 3 通道'),('#187896','高 q 合并 5 通道')]]
        if field=='low_count_apex_eV':handles.append(Line2D([],[],marker='.',ls='none',color='#bbbbbb',label='总谱峰顶（300 帧）'))
        fig.legend(handles=handles,loc='outside lower center',ncol=4)
        fig.suptitle(title+'：逐谱双成分拟合，无跨 q 上升约束\n'
            '实心：条件模型内稳定且留出预测改善；空心：候选；叉号：未分辨；竖线：95% 帧重采样区间',fontsize=13)
        save(fig,output,name)
    tasklookup={t['key']:t for t in tasks}
    for sid in ['590','n0','20w']:
        fig,axes=plt.subplots(1,5,figsize=(17,4.5),layout='constrained')
        group=rr(sid,5,1)
        for ax,target in zip(axes,[.001,.0025,.005,.010,.015]):
            r=min(group,key=lambda r:abs(r['q_Ainv']-target));f=fits[r['key']];a=f['A'];s=f['S']
            d=local_data(tasklookup[r['key']]);e=d['energy']/1000;band=(e>=.3)&(e<=1.8)
            test=d['frames'][0,250:].astype(float).sum(0)/a['test_scale']
            sigma=d['sigma_test'][0]/a['test_scale'];comp=np.asarray(a['components'])[0]
            norm=max(np.asarray(a['prediction'])[0,band].max(),1e-100)
            ax.errorbar(e[band],test[band]/norm,yerr=sigma[band]/norm,fmt='.',ms=2.5,
                        color='#777777',alpha=.6,elinewidth=.5,label='留出 50 帧（剂量折算）')
            ax.plot(e[band],np.asarray(a['prediction'])[0,band]/norm,color='#1c6e91',lw=1.6,label='双成分总谱')
            ax.plot(e[band],np.asarray(s['prediction'])[0,band]/norm,color='#b3b3b3',ls='--',lw=1,label='单 B1 成分总谱')
            for j,color,label in [(0,'#be7535','低 E₀ 成分'),(1,'#258b70','高 E₀ 成分')]:
                ax.plot(e[band],comp[band,j]/norm,color=color,lw=1.6,label=label)
            ax.plot(e[band],comp[band,2:].sum(1)/norm,':',color='#8f7298',lw=1.2,label='邻近带与背景')
            ax.axvline(r['low_count_apex_eV'],color='#be7535',ls=':',lw=.8)
            ax.set(title=f'q={r["q_Ainv"]:.4f} | {r["count_status"]}\n'
                f'低峰顶 {r["low_count_apex_eV"]:.3f}；E₀ {r["low_E0_eV"]:.3f} eV',
                xlabel='能量损失 / eV',ylabel='计数 / 拟合总谱最大值',ylim=(0,1.25),xlim=(.3,1.8))
        handles,labels=axes[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=9)
        fig.suptitle(f'{sid}：低能成分单独展示；前 250 帧拟合，后 50 帧检验\n'
                     '实际接收积分、运动学权重与能量卷积；每谱独立；图中点未参与参数拟合',fontsize=13)
        save(fig,output,f'03_components_{sid}')
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
    for ax,sid in zip(axes,['590','n0','20w']):
        for field,color,label in [('low_count_apex_eV','#be7535','低能成分计数峰顶'),('low_E0_eV','#1c6e91','低能振子 E₀')]:
            group=[r for r in rows if r['session']==sid and r['high_channels']==5 and r['q_Ainv']>0]
            anchor=[c for c in checks if c['session']==sid and c['q_Ainv']>0 and c['variant']=='window400']
            for c in anchor:
                r=next(r for r in group if r['key']==c['key']);q=r['q_Ainv']
                ax.plot([q,q],[r[field],c[field]],color=color,lw=1.1,alpha=.6)
                ax.plot(q,r[field],'o',ms=4,color=color)
                ax.plot(q,c[field],'s',ms=4,color=color,mfc='white')
            ax.plot([],[],color=color,label=label)
        ax.set(title=sid,xlabel='正 q / Å⁻¹',ylabel='能量 / eV',ylim=(.3,1.6));ax.grid(alpha=.15)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=2)
    fig.suptitle('能窗敏感性：圆点从 300 meV 起拟合，方点从 400 meV 起拟合\n'
                 '其余数据与模型相同；高 q 的大幅分歧是分峰不唯一的直接信号',fontsize=13)
    save(fig,output,'04_window_sensitivity')
    dump(output/'figures_provenance.json',dict(code_sha256=digest(Path(__file__)),
         rows_sha256=digest(output/'low_component_points.csv'),figures={f.name:digest(f) for f in output.glob('0*.png')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--ridge',type=Path,required=True);a=p.parse_args();summarize(a.output,a.ridge)
