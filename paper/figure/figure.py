import HamiltonianPy as HP
import numpy as np
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import pdb

def lattice():
    plt.ion()
    gs=plt.GridSpec(21,2)
    gs.update(left=0.02,right=0.98,top=1.00,bottom=0.09,hspace=0.25,wspace=0.20)

    ax=plt.subplot(gs[0:7,0])
    ax.axis('off')
    ax.set_xlim(-0.1,5.6)
    ax.set_ylim(-1.0,2.25)
    PA,PB,t=np.array([0.0,0.0]),np.array([0.5,0.8]),np.array([1.0,0.0])
    lattice=HP.Lattice(name='FBFM',rcoords=HP.tiling(cluster=[PA,PB],vectors=[t],translations=xrange(6)),neighbours=2)
    for bond in lattice.bonds:
        p1,p2=bond.spoint,bond.epoint
        if bond.neighbour==0:
            if p1.pid.site<len(lattice)-1:
                alpha=1.0 if p1.pid.site in (0,1,2,3,4,5,6,9,10) else 0.5
                ax.scatter(p1.rcoord[0],p1.rcoord[1],edgecolors='none',color='blue' if p1.pid.site%2==0 else 'red',alpha=alpha,zorder=2)
        else:
            if not (bond.neighbour==2 and p1.pid.site%2==1 and p2.pid.site%2==1 or p1.pid.site==len(lattice)-1 or p2.pid.site==len(lattice)-1):
                ax.plot([p1.rcoord[0],p2.rcoord[0]],[p1.rcoord[1],p2.rcoord[1]],lw=2,color='purple' if bond.neighbour==1 else 'green',zorder=1)
    ax.text(0.90,-0.70,'A',fontsize=16,color='black')
    ax.text(1.40,+1.00,'B',fontsize=16,color='black')
    ax.annotate(s='',xy=[2.0,0.0],xytext=[3.0,0.0],arrowprops={'color':'green','linewidth':2,'arrowstyle':'<->','zorder':3})
    ax.annotate(s='',xy=[2.0,0.0],xytext=[2.5,0.8],arrowprops={'color':'purple','linewidth':2,'arrowstyle':'<->','zorder':3})
    ax.annotate(s='',xy=[2.5,0.8],xytext=[3.0,0.0],arrowprops={'color':'purple','linewidth':2,'arrowstyle':'<->','zorder':3})
    ax.text(2.41,-0.62,'$t$',fontsize=18,color='black')
    ax.text(2.00,0.40,'$\lambda$',fontsize=18,color='black')
    ax.text(2.75,0.40,'$\lambda$',fontsize=18,color='black')
    ax.text(2.39,0.95,'$\epsilon$',fontsize=18,color='black')
    dy=0.8
    ax.annotate(s='',xy=[3.45,0.8-dy/2],xytext=[3.45,0.8+dy/2+0.1],arrowprops={'color':'red','linewidth':2,'arrowstyle':'<-','zorder':3})
    ax.annotate(s='',xy=[3.55,0.8-dy/2-0.1],xytext=[3.55,0.8+dy/2],arrowprops={'color':'red','linewidth':2,'arrowstyle':'->','zorder':3})
    ax.annotate(s='',xy=[3.95,-dy/2],xytext=[3.95,dy/2+0.1],arrowprops={'color':'blue','linewidth':2,'arrowstyle':'<-','zorder':3})
    ax.annotate(s='',xy=[4.05,-dy/2-0.1],xytext=[4.05,dy/2],arrowprops={'color':'blue','linewidth':2,'arrowstyle':'->','zorder':3})
    ax.text(3.35,+1.24,'$U_d$',fontsize=18,color='black')
    ax.text(3.90,-0.97,'$U_s$',fontsize=18,color='black')
    ax.text(-0.1,1.47,"(a)",ha='left',fontsize=18,color='black')

    ax=plt.subplot(gs[7:,0])
    ax.axis('off')
    ax.set_xlim(-1.5,6.1)
    ax.set_ylim(-0.4,4.5)
    sp,ep,t=np.array([0.0]),np.array([0.75]),np.array([1.0])
    lattice=HP.Lattice(name='FBFM',rcoords=HP.tiling(cluster=[sp,ep],vectors=[t],translations=xrange(6)),neighbours=2)
    for bond in lattice.bonds:
        if bond.neighbour==2:
            x1,x2,dy=bond.spoint.rcoord[0],bond.epoint.rcoord[0],0.8
            ax.plot([x1,x2],[0.0,0.0],lw=2,color='black')
            ax.plot([x1,x2],[3.0,3.0],lw=2,color='black')
            ax.annotate(s='',xy=[(x1+x2)/2,-dy/2],xytext=[(x1+x2)/2,dy/2+0.1],arrowprops={'color':'green','linewidth':2,'arrowstyle':'<-','zorder':3})
            if np.abs((x1+x2)/2-4.375)<HP.RZERO:
                ax.annotate(s='',xy=[(x1+x2)/2-0.05,3.0-dy/2],xytext=[(x1+x2)/2-0.05,3.0+dy/2+0.1],arrowprops={'color':'green','linewidth':2,'arrowstyle':'<-','zorder':3})
                ax.annotate(s='',xy=[(x1+x2)/2+0.05,3.0-dy/2-0.1],xytext=[(x1+x2)/2+0.05,3.0+dy/2],arrowprops={'color':'green','linewidth':2,'arrowstyle':'->','zorder':3})
            elif np.abs((x1+x2)/2-1.375)>HP.RZERO:
                ax.annotate(s='',xy=[(x1+x2)/2,3.0-dy/2],xytext=[(x1+x2)/2,3.0+dy/2+0.1],arrowprops={'color':'green','linewidth':2,'arrowstyle':'<-','zorder':3})
    ax.text(1.25,-0.9,'$k_i$',fontsize=18,color='black')
    ax.text(4.25,-0.9,'$k_j$',fontsize=18,color='black')
    ax.text(1.25,2.1,'$k_i$',fontsize=18,color='black')
    ax.text(4.25,2.1,'$k_j$',fontsize=18,color='black')
    ax.plot([1.375,1.375],[3.5,4.3],ls='--',color='black')
    ax.plot([4.375,4.375],[3.5,4.3],ls='--',color='black')
    ax.annotate('',xy=[1.375,3.7],xytext=[4.375,3.7],arrowprops={'color':'black','linewidth':1.5,'arrowstyle':'<->','zorder':3})
    ax.text(2.875,3.7,'$q=k_i-k_j$',ha='center',va='bottom',fontsize=18,color='black')
    ax.annotate('',xy=[2.875,0.5],xytext=[2.875,2.5],arrowprops={'color':'black','linewidth':1.5,'arrowstyle':'<-','zorder':3})
    ax.text(2.975,1.5,'$\Delta S_z=-1$',ha='left',va='center',fontsize=18,color='black')
    ax.text(-0.1,0.0,'$\left|\mathrm{GS}\\rightangle\\right.$:',ha='right',va='center',fontsize=14,color='black')
    ax.text(-0.1,3.0,'$S_q^-\left|\mathrm{GS}\\rightangle\\right.$:',ha='right',va='center',fontsize=14,color='black')
    ax.text(-1.5,4.3,'(c)',ha='left',va='center',fontsize=18,color='black')

    ax=plt.subplot(gs[1:20,1])
    result=np.loadtxt('../../result/tba/1DIF_S2x(1P-1O)_A_1.0_1.4_-0.04_TBA_EB.dat')
    ax.plot(result[:,0],result[:,1:],color='green',lw=2)
    ax.minorticks_on()
    ax.set_ylim(-2.5,4.5)
    ax.set_xticks(np.linspace(0,400,5))
    ax.set_xticklabels(['-$\pi$','','$0$','','$\pi$'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel("k",fontdict={'fontsize':18})
    ax.set_ylabel("$E/t$",fontdict={'fontsize':18})
    ax.text(-80,4.0,"(b)",fontsize=18,color='black')

    pdb.set_trace()
    plt.savefig('lattice.pdf')
    plt.close()

def spectrum():
    plt.ion()
    fig,axes=plt.subplots(nrows=1,ncols=3)
    fig.subplots_adjust(left=0.10,right=0.98,top=0.96,bottom=0.14,hspace=0.2,wspace=0.25)

    start,inc,vmin,vmax=11,-3,0,12
    cmap=cmx.ScalarMappable(norm=colors.Normalize(vmin=vmin,vmax=vmax),cmap=plt.get_cmap('viridis'))
    for i,parameter in enumerate([0.02,0.409,0.78]):
        names=[
            '../../result/ed/1DIF_S2x(8P-1O)_FSTR(32,8,3.0)_1.0_1.4_-0.04_%s_1.0_TrFED_TREDEB.dat'%parameter,
            '../../result/fbfm/1DIF_S2x(1P-1O)_up_1.0_1.4_-0.04_%s_1.0_FBFM_EB8.dat'%parameter,
            '../../result/fbfm/1DIF_S2x(1P-1O)_up_1.0_1.4_-0.04_%s_1.0_FBFM_EB16.dat'%parameter,
            '../../result/fbfm/1DIF_S2x(1P-1O)_up_1.0_1.4_-0.04_%s_1.0_FBFM_EB24.dat'%parameter,
            '../../result/fbfm/1DIF_S2x(1P-1O)_up_1.0_1.4_-0.04_%s_1.0_FBFM_EB60.dat'%parameter,
            '../../result/fbfm/1DIF_S2x(1P-1O)_up_1.0_1.4_-0.04_%s_1.0_FBFM_EB800.dat'%parameter,
        ]
        for j,(name,label) in enumerate(zip(names,('','$N_q=8$','$N_q=16$','$N_q=24$','$N_q=60$','$N_q=800$'))):
            result=np.loadtxt(name)
            if j==0:
                result=result[np.array([4,5,6,7,0,1,2,3,4]),:]
                xs=np.array(range(result.shape[0]))*1.0/(result.shape[0]-1)
                axes[i].plot(xs,(result[:,1:]+16.0)/parameter,'.',color='red',lw=8,zorder=2,clip_on=False)
            else:
                result=result[np.array(range(result.shape[0])+[0]),:]
                xs=np.array(range(result.shape[0]))*1.0/(result.shape[0]-1)
                if i==1:
                    axes[i].plot(xs,result[:,1]/parameter,color=cmap.to_rgba(start+(j-1)*inc),ls='-',lw=1.5,zorder=1)
                    axes[i].plot(xs,result[:,2]/parameter,color=cmap.to_rgba(start+(j-1)*inc),ls='-',label=label,lw=1.5,zorder=1)
                else:
                    axes[i].plot(xs,result[:,1:]/parameter,color=cmap.to_rgba(start+(j-1)*inc),ls='-',label=label,lw=1.5,zorder=1)
            if i==0 and label=='$N_q=800$':
                from mpl_toolkits.axes_grid.inset_locator import inset_axes
                ax=inset_axes(axes[0],width="40%",height=1.0,loc=1)
                ax.plot(xs,result[:,1:]/parameter,ls='-',color=cmap.to_rgba(start+(j-1)*inc),lw=1.5,zorder=1)
                ax.minorticks_on()
                ax.set_xlim(0.25,0.75)
                ax.set_ylim(0.5,0.6)
                ax.set_xticks([0.25,0.5,0.75])
                ax.set_xticklabels(['','$\pi$',''])
                ax.set_yticks([0.52,0.58])
                for tick in ax.get_xticklabels():
                    tick.set_fontsize(14)
                for tick in ax.get_yticklabels():
                    tick.set_fontsize(14)
        axes[i].minorticks_on()
        axes[i].set_xlim(-0.01,1.01)
        axes[i].set_xticks(np.linspace(0.0,1.0,5))
        axes[i].set_xticklabels(['0','','$\pi$','','$2\pi$'])
        axes[i].set_xlabel('q',fontdict={'fontsize':18})
        if i==0:
            axes[i].set_ylim(-0.01,2.0)
            axes[i].set_yticks(np.linspace(0.0,2.0,5))
            axes[i].set_ylabel('$E/U_s$',fontdict={'fontsize':18})
            axes[i].text(0.1,1.8,"(a)",fontsize=18,color='black')
        if i==1:
            axes[i].set_ylim(-0.003,0.6)
            axes[i].set_yticks(np.linspace(0.0,0.6,3))
            axes[i].text(0.1,0.54,"(b)",fontsize=18,color='black')
            leg=axes[i].legend(loc='lower center',fancybox=True,shadow=False,prop={'size': 14})
            leg.get_frame().set_alpha(0.5)
        if i==2:
            axes[i].set_ylim(-0.0025,0.5)
            axes[i].set_yticks(np.linspace(0.0,0.5,6))
            axes[i].text(0.1,0.45,"(c)",fontsize=18,color='black')
        for tick in axes[i].get_xticklabels():
            tick.set_fontsize(18)
        for tick in axes[i].get_yticklabels():
            tick.set_fontsize(18)

    pdb.set_trace()
    plt.savefig('spectrum.pdf')
    plt.close()

def gap_berry_and_phase_diagram():
    plt.ion()
    fig,axes=plt.subplots(nrows=1,ncols=2)
    fig.subplots_adjust(left=0.10,right=0.98,top=0.96,bottom=0.17,hspace=0.2,wspace=0.40)

    ax=axes[0]
    twinx=ax.twinx()
    fbbp=np.loadtxt('../../result/fbfm/S2x_BP_FBFM.dat')
    twinx.plot(fbbp[:,0],np.abs(fbbp[:,1]),'-',color='green',lw=2,zorder=1)
    twinx.plot(fbbp[:,0],np.abs(fbbp[:,2]),'--',color='green',lw=2,zorder=1)
    twinx.minorticks_on()
    twinx.set_ylim(-0.001,1.15)
    twinx.set_yticks(np.linspace(0,1.0,3))
    twinx.set_yticklabels(['0','$\\frac{\pi}{2}$','$\pi$'])
    for tick in twinx.get_yticklabels():
        tick.set_fontsize(18)
        tick.set_color('green')
    twinx.set_ylabel('Berry Phase',color='green',fontdict={'fontsize':18})

    fbgap=np.loadtxt('../../result/fbfm/S2x_GAP_FBFM_EX.dat')
    ax.plot(fbgap[:,0],fbgap[:,-1],color='blue',ls='-',lw=2,zorder=1)
    ax.minorticks_on()
    ax.set_xlim(0.0,0.8)
    ax.set_ylim(0.0,0.1)
    ax.set_xticks(np.linspace(0.0,0.8,5))
    ax.set_xticklabels(['0.0','0.2','0.4','0.6','0.8'])
    ax.set_yticks(np.linspace(0.0,0.1,6))
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
        tick.set_color('blue')
    ax.set_xlabel("$U_s/U_d$",fontdict={'fontsize':22})
    ax.set_ylabel("$\Delta_\pi$",color='blue',fontdict={'fontsize':20})
    ax.text(0.05,0.09,"(a)",fontsize=18,color='black',zorder=4)

    ax=axes[1]
    xs=np.array([0.02,0.1,0.2,0.3,0.409,0.5,0.6,0.7,0.8,0.9,1.0])
    ys=np.array([1.80,1.68,1.56,1.47,1.40,1.34,1.29,1.25,1.22,1.19,1.16])
    ax.plot(xs,ys,lw=2)
    ax.plot([0.02,0.78],[1.4,1.4],ls='--',color='grey',dashes=(0.5,5.0),dash_capstyle='round',lw=2)
    ax.minorticks_on()
    ax.set_xlim(0.0,0.8)
    ax.set_ylim(1.0,1.8)
    ax.set_xticks(np.linspace(0.0,0.8,5))
    ax.set_xticklabels(['0.0','0.2','0.4','0.6','0.8'])
    ax.set_yticks(np.linspace(1.0,1.8,5))
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel("$U_s/U_d$",fontdict={'fontsize':22})
    ax.set_ylabel("$\lambda/t$",fontdict={'fontsize':20})
    ax.text(0.20,1.20,"FM",fontsize=22,color='black')
    ax.text(0.40,1.50,"TFM",fontsize=22,color='black')
    ax.text(0.7,1.70,"(b)",fontsize=18,color='black',zorder=4)

    pdb.set_trace()
    plt.savefig('gap_berry_and_phase_diagram.pdf')
    plt.close()

def edge():
    plt.ion()
    gs=plt.GridSpec(3,3)
    gs.update(left=0.14,right=0.96,top=0.97,bottom=0.11,hspace=0.55,wspace=0.20)

    ax=plt.subplot(gs[0:2,:])
    result=np.loadtxt('../../result/fbfm/1DIF_S2x(60O-1O)_up_1.0_1.4_-0.04_1.0_FBFM_EDGE.dat')
    xs=result[:,0]
    ax.plot(xs,result[:,1:59]/xs[:,np.newaxis],color='grey',lw=1,alpha=0.5,zorder=1)
    ax.plot(xs,result[:,59]/xs,color='blue',lw=2,zorder=2)
    ax.plot(xs,result[:,60]/xs,color='green',lw=2,zorder=2)
    ax.plot(xs,result[:,61]/xs,color='purple',lw=2,zorder=2)
    ax.plot(xs,result[:,62]/xs,color='red',lw=2,zorder=2)
    ax.plot(xs,result[:,63:110]/xs[:,np.newaxis],color='grey',lw=2,alpha=0.5,zorder=1)

    ax.minorticks_on()
    ax.set_ylim(0.19,0.6)
    ax.set_xticks(np.linspace(0.0,0.8,5))
    ax.set_yticks(np.linspace(0.2,0.6,5))
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel("$U_s/U_d$",fontdict={'fontsize':20})
    ax.set_ylabel("$E/Us$",fontdict={'fontsize':20})
    ax.text(0.7,0.56,"(a)",fontsize=18,color='black')

    for i,(parameter,tag) in enumerate(zip(['0.02','0.3','0.78'],['(b)','(c)','(d)'])):
        ax=plt.subplot(gs[2,i])
        result=np.loadtxt('../../result/fbfm/1DIF_S2x(60O-1O)_up_1.0_1.4_-0.04_%s_1.0_FBFM_POS.dat'%parameter)
        ax.plot(result[:,0],result[:,2].real,color='blue',lw=2)
        ax.plot(result[:,0],result[:,3].real,color='green',lw=2)
        ax.plot(result[:,0],result[:,4].real,color='purple',lw=2)
        ax.plot(result[:,0],result[:,5].real,color='red',lw=2)
        ax.minorticks_on()
        ax.set_ylim(0.0,0.45)
        ax.set_xticks(np.linspace(0,120,4))
        ax.set_yticks(np.linspace(0.0,0.4,3))
        for tick in ax.get_xticklabels():
            tick.set_fontsize(16)
        ax.set_xlabel('position',fontdict={'fontsize':20})
        if i==0:
            for tick in ax.get_yticklabels():
                tick.set_fontsize(18)
            ax.set_ylabel('$\Delta\langle S_z\\rangle$',fontdict={'fontsize':18})
        else:
            ax.set_yticklabels(['']*3)
        ax.text(10,0.32,tag,fontsize=18,color='black')

    pdb.set_trace()
    plt.savefig('edge.pdf')
    plt.close()

def delta_spectrum():
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    plt.ion()
    gs=plt.GridSpec(7,2)
    gs.update(left=0.13,right=0.97,top=0.99,bottom=0.07,hspace=0.9,wspace=0.1)

    for i,(parameter,tag) in enumerate(zip([-0.04,0.18,0.26,0.41],['(a)','(b)','(c)','(d)'])):
        data=np.load('../../result/fbfm/1DIF_S2x(1P-1O)_up_1.0_1.4_%s_0.3_1.0_FBFM_EB800.npz'%HP.decimaltostr(parameter))['data']
        ax=plt.subplot(gs[0:2,i%2]) if i<2 else plt.subplot(gs[2:4,i%2])
        ax.plot(data[:,0],data[:,1:3],color='green',lw=2,zorder=2)
        ax.plot(data[:,0],data[:,3:],color='grey',lw=2,alpha=0.3,zorder=1)

        ax.minorticks_on()
        ax.set_ylim(0.0,0.25)
        ax.set_xticks(np.linspace(0,800,5))
        ax.set_xticklabels(['0','','$\pi$','','$2\pi$'])
        ax.set_xlabel('q',fontdict={'fontsize':18})
        ax.set_yticks(np.linspace(0.0,0.25,6))
        ax.set_yticklabels(['0','','0.1','','0.2',''] if i in (0,2) else ['']*6)
        if i in (0,2): ax.set_ylabel('$E/U_s$',fontdict={'fontsize':18})
        for tick in ax.get_xticklabels():
            tick.set_fontsize(18)
        for tick in ax.get_yticklabels():
            tick.set_fontsize(18)
        ax.text(50,0.21,tag,fontsize=22,color='black')

        ax=inset_axes(ax,width="30%",height=1.0,loc=1)
        data=np.loadtxt('../../result/tba/1DIF_S2x(1P-1O)_A_1.0_1.4_%s_TBA_EB.dat'%HP.decimaltostr(parameter))
        ax.plot(data[:,0],data[:,1:],lw=1.5,color='blue',zorder=1)
        ax.minorticks_on()
        ax.set_ylim(-2.5,4.5)
        ax.set_xticks([0,200,400])
        ax.set_xticklabels(['$-\pi$','k','$\pi$'])
        ax.set_yticks([-2.5,1.0,4.5])
        ax.set_yticklabels(['-2.5','$E/t$','4.5'])
        for tick in ax.get_xticklabels():
            tick.set_fontsize(14)
            tick.set_color('blue')
        for tick in ax.get_yticklabels():
            tick.set_fontsize(14)
            tick.set_color('blue')

    ax=plt.subplot(gs[4:,:])
    xs1=np.array([0.02,0.10,0.15,0.20,0.250,0.30,0.350,0.360,0.370,0.380,0.39,0.395,0.40,0.405,0.41])
    ys1=np.array([0.36,0.33,0.31,0.28,0.255,0.22,0.175,0.165,0.155,0.145,0.13,0.120,0.11,0.100,0.00])
    xs2=np.array([0.02,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50,0.60,0.70,0.80])
    ys2=np.array([0.38,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51])
    ax.plot(xs1,ys1,lw=2)
    ax.plot(xs2,ys2,lw=2)
    ax.scatter([0.3,0.3,0.3,0.3],[0.0,0.22,0.30,0.45],s=100,color='red',marker='*',edgecolors='none',zorder=4,clip_on=False,alpha=0.8)
    ax.minorticks_on()
    ax.set_xlim(0.1,0.7)
    ax.set_ylim(0.0,0.6)
    ax.set_xticks(np.linspace(0.1,0.7,4))
    ax.set_xticklabels(['0.1','0.3','0.5','0.7'])
    ax.set_yticks(np.linspace(0.0,0.6,4))
    ax.set_yticklabels(['0','0.2','0.4','0.6'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel("$U_s/U_d$",fontdict={'fontsize':20})
    ax.set_ylabel("$\Delta/t$",fontdict={'fontsize':20})
    ax.text(0.20,0.10,"FM",fontsize=22,color='black')
    ax.text(0.45,0.20,"TFM",fontsize=22,color='black')
    ax.text(0.15,0.48,"NFM",fontsize=22,color='black')
    ax.text(0.60,0.52,"(e)",fontsize=22,color='black',zorder=4)

    pdb.set_trace()
    plt.savefig('delta_spectrum.pdf')
    plt.close()

def v_spectrum():
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    plt.ion()
    gs=plt.GridSpec(7,2)
    gs.update(left=0.13,right=0.97,top=0.98,bottom=0.07,hspace=0.9,wspace=0.1)

    for i,(parameter,tag) in enumerate(zip([0.0,0.09,0.4,0.64],['(a)','(b)','(c)','(d)'])):
        data=np.load('../../result/fbfm/1DIF_1.0_1.25_-0.4375_0.3_1.0_%s_EB800.npz'%HP.decimaltostr(parameter))['data']
        ax=plt.subplot(gs[0:2,i%2]) if i<2 else plt.subplot(gs[2:4,i%2])
        ax.plot(data[:,0],data[:,1:3],color='green',lw=2,zorder=2)
        ax.plot(data[:,0],data[:,3:],color='grey',lw=2,alpha=0.3,zorder=1)

        ax.minorticks_on()
        ax.set_ylim(0.0,0.30)
        ax.set_xticks(np.linspace(0,800,5))
        ax.set_xticklabels(['0','','$\pi$','','$2\pi$'])
        ax.set_xlabel('q',fontdict={'fontsize':18})
        ax.set_yticks(np.linspace(0.0,0.3,4))
        ax.set_yticklabels(['0','0.1','0.2','0.3'] if i in (0,2) else ['']*6)
        if i in (0,2): ax.set_ylabel('$E/U_s$',fontdict={'fontsize':18})
        for tick in ax.get_xticklabels():
            tick.set_fontsize(18)
        for tick in ax.get_yticklabels():
            tick.set_fontsize(18)
        ax.text(50,0.25,tag,fontsize=22,color='black')

    ax=plt.subplot(gs[4:,:])
    xs1=np.array([0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.71])
    ys1=np.array([0.05,0.08,0.09,0.09,0.08,0.05,0.01,0.00])
    xs2=np.array([0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90])
    ys2=np.array([0.38,0.54,0.64,0.72,0.81,0.85,0.89,0.93,0.97])
    ax.plot(xs1,ys1,lw=2)
    ax.plot(xs2,ys2,lw=2)
    ax.scatter([0.3,0.3,0.3,0.3],[0.00,0.09,0.40,0.64],s=100,color='red',marker='*',edgecolors='none',zorder=4,clip_on=False,alpha=0.8)
    ax.minorticks_on()
    ax.set_xlim(0.1,0.9)
    ax.set_ylim(0.0,1.0)
    ax.set_xticks(np.linspace(0.1,0.9,5))
    ax.set_xticklabels(['0.1','0.3','0.5','0.7','0.9'])
    ax.set_yticks(np.linspace(0.0,1.0,6))
    ax.set_yticklabels(['0','0.2','0.4','0.6','0.8','1.0'])
    for tick in ax.get_xticklabels():
        tick.set_fontsize(18)
    for tick in ax.get_yticklabels():
        tick.set_fontsize(18)
    ax.set_xlabel("$U_s/U_d$",fontdict={'fontsize':20})
    ax.set_ylabel("$V/U_d$",fontdict={'fontsize':20})
    ax.text(0.35,0.01,"FM",fontsize=22,color='black')
    ax.text(0.50,0.35,"TFM",fontsize=22,color='black')
    ax.text(0.15,0.75,"NFM",fontsize=22,color='black')
    ax.text(0.80,0.83,"(e)",fontsize=22,color='black',zorder=4)

    pdb.set_trace()
    plt.savefig('v_spectrum.pdf')
    plt.close()

if __name__=='__main__':
    for arg in sys.argv:
        if arg in ('1','all'): lattice()
        if arg in ('2','all'): spectrum()
        if arg in ('3','all'): gap_berry_and_phase_diagram()
        if arg in ('4','all'): edge()
        if arg in ('5','all'): delta_spectrum()
        if arg in ('6','all'): v_spectrum()
