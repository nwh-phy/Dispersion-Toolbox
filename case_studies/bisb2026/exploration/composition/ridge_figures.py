"""Display observed adaptive-q apices, including unstable extractions."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from ridge_extract import smooth_counts
from kinematic_audit import save
from io_utils import digest,dump


def make(output):
    font=Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
    if font.exists():
        font_manager.fontManager.addfont(font)
        plt.rcParams['font.family']=font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    rows=list(csv.DictReader((output/'peaks.csv').open()))
    data={}
    for sid in ['590','n0','20w']:
        for n in [3,5]:
            manifest=json.loads((output/f'bin{n}/manifest.json').read_text())
            expected=next(x for x in manifest['sessions'] if x['session']==sid)
            assert digest(output/f'bin{n}/{sid}.npz')==expected['prepared_sha256']
            data[sid,n]=dict(np.load(output/f'bin{n}/{sid}.npz'))

    def get_rows(sid,n,sign=None):
        return [r for r in rows if r['session']==sid and int(r['high_channels'])==n
                and (sign is None or float(r['q_Ainv'])*sign>0)]

    for sessions,name in [(['590','n0'],'01_adaptive_ridge_10w'),(['20w'],'02_adaptive_ridge_20w')]:
        fig,axes=plt.subplots(len(sessions),2,figsize=(13,4.2*len(sessions)+.6),
                             layout='constrained',squeeze=False)
        for row,sid in zip(axes,sessions):
            for ax,n in zip(row,[3,5]):
                d=data[sid,n];e=d['energy']/1000;band=(e>=.3)&(e<=1.8)
                y=d['frames'].astype(float).sum(1)[:,band]
                y/=y.max(1)[:,None];y[abs(d['q'])<1e-10]=np.nan
                edges=np.r_[d['q']-d['width']/2,d['q'][-1]+d['width'][-1]/2]
                ee=np.r_[e[band]-.008,e[band][-1]+.008]
                im=ax.pcolormesh(edges,ee,y.T,cmap='turbo',vmin=.3,vmax=1,rasterized=True)
                for r in get_rows(sid,n):
                    q=float(r['q_Ainv']);p=float(r['peak_meV'])/1000
                    if not np.isfinite(p):continue
                    lo=float(r['ci_lo_meV'])/1000;hi=float(r['ci_hi_meV'])/1000
                    ax.plot([q,q],[lo,hi],color='black',alpha=.55,lw=.65)
                    ax.plot(q,p,'o' if r['stable']=='True' else 'x',ms=4,
                            color='white',mec='#202020' if r['stable']=='True' else 'white',mew=.7)
                for s in [-1,1]:ax.axvline(s*(.006+float(d['dq'])/2),ls='--',color='white',lw=.8,alpha=.8)
                ax.set(title=f'{sid} | 小 q 单通道；高 q 合并 {n} 通道',xlabel='q / Å⁻¹',
                       ylabel='能量损失 / eV',ylim=(.3,1.8),xlim=(-.02,.02))
        fig.colorbar(im,ax=axes,shrink=.8,label='真实计数 / 每列最大值')
        fig.suptitle('B1 实测峰顶随 |q| 上升：逐谱独立提取，没有施加单调上升约束\n'
                     '热图未作能量平滑、未扣背景；○ 较稳定，× 需复核；竖线为 95% 帧重采样区间',fontsize=13)
        save(fig,output,name)

    fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    for col,sid in enumerate(['590','n0','20w']):
        for row,sign in enumerate([1,-1]):
            ax=axes[row,col]
            for n,color in [(3,'#c27839'),(5,'#176f93')]:
                rr=sorted(get_rows(sid,n,sign),key=lambda r:abs(float(r['q_Ainv'])))
                for r in rr:
                    q=abs(float(r['q_Ainv']));p=float(r['peak_meV'])/1000
                    if not np.isfinite(p) or (n==5 and q<=.006+1e-10):continue
                    native=q<=.006+1e-10;c='#333333' if native else color
                    lo=float(r['ci_lo_meV'])/1000;hi=float(r['ci_hi_meV'])/1000
                    ax.plot([q,q],[lo,hi],color=c,lw=.7,alpha=.6)
                    ax.plot([q-float(r['width_Ainv'])/2,q+float(r['width_Ainv'])/2],[p,p],color=c,lw=.6)
                    ax.plot(q,p,'o',ms=3.8,color=c,mfc=c if r['stable']=='True' else 'white',mew=.9)
            ax.axvspan(0,.006,color='#888888',alpha=.08)
            ax.set(title=f'{sid} | {"正 q" if sign>0 else "负 q"}',xlabel='|q| / Å⁻¹',
                   ylabel='实测峰顶 / eV',xlim=(0,.02),ylim=(.3,1.85))
            ax.grid(alpha=.15)
    handles=[Line2D([],[],marker='o',color=c,ls='none',label=l) for c,l in
             [('#333333','小 q：单通道'),('#c27839','高 q：3 通道'),('#176f93','高 q：5 通道')]]
    fig.legend(handles=handles,loc='outside lower center',ncol=3)
    fig.suptitle('上升幅度约 1 eV；高 q 转为宽平台，部分峰位随噪声跳动\n'
                 '空心点需复核；横线表示 q 积分宽度，竖线为 95% 重采样区间；逐点结果未经曲线约束',fontsize=13)
    save(fig,output,'03_signed_peak_positions')

    fig,axes=plt.subplots(2,5,figsize=(16,7),layout='constrained')
    for sid,axs in zip(['590','n0'],axes):
        d=data[sid,5];e=d['energy']/1000
        for t,ax in zip([.001,.0025,.005,.010,.015],axs):
            qi=abs(d['q']-t).argmin();q=d['q'][qi];y=d['frames'][qi].astype(float).sum(0)
            s=smooth_counts(d['energy'],y,100);band=(e>=.3)&(e<=1.8);scale=s[band].max()
            rr=min(get_rows(sid,5),key=lambda r:abs(float(r['q_Ainv'])-q));p=float(rr['peak_meV'])/1000
            ax.plot(e[band],y[band]/scale,lw=.8,color='#adbbc4',label='实测 16 meV 通道')
            ax.plot(e[band],s[band]/scale,lw=1.8,color='#176f93',label='100 meV FWHM 平滑')
            ax.axvline(p,color='#b75c3b',ls='--',lw=1,label='独立提取峰顶')
            ax.set(title=f'{sid} | q={q:.4f}\nN={rr["native_channels"]}；峰顶 {p:.3f} eV',
                   xlabel='能量损失 / eV',ylabel='归一化计数',ylim=(0,1.15),xlim=(.3,1.8))
            ax.grid(alpha=.12)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=3)
    fig.suptitle('原始单谱支持峰顶上升：从约 0.45 eV 到 1.45–1.48 eV\n'
                 '300 帧漂移校正后求和；没有扣背景、除运动学因子或设定两个振子',fontsize=13)
    save(fig,output,'04_real_spectra_with_peaks')
    dump(output/'figures_provenance.json',dict(code_sha256=digest(Path(__file__)),
         peaks_sha256=digest(output/'peaks.csv'),figures={p.name:digest(p) for p in output.glob('0*.png')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    make(p.parse_args().output)
