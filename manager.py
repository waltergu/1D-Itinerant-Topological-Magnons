from HamiltonianPy import *
from source import *
from collections import OrderedDict
from fractions import Fraction
import numpy as np
import mkl

def tbatasks(parameters,lattice,terms,tbaidfmap=tbaidfa,jobs=('EB',)):
    import HamiltonianPy.FreeSystem as TBA
    tba=tbaconstruct(parameters,lattice,terms,tbaidfmap)
    if 'EB' in jobs:
        path=KSpace(reciprocals=lattice.reciprocals,segments=[(-0.5,0.5)],end=True,nk=401) if len(lattice.vectors)>0 else None
        tba.register(EB(name='EB',path=path,run=TBA.TBAEB))
        tba.summary()
    if 'EDGE' in jobs:
        from scipy.linalg import eigh
        import matplotlib.pyplot as plt
        import pdb
        m=len(lattice)
        es,evs=eigh(tba.matrix(),eigvals_only=False)
        plt.ion()
        fig,axes=plt.subplots(nrows=2,ncols=1)
        fig.subplots_adjust(left=0.10,right=0.98,top=0.96,bottom=0.14,hspace=0.2,wspace=0.25)
        axes[0].plot([0,1],np.vstack([es[:m/2-1],es[:m/2-1]]),color='black',lw=1,zorder=1)
        axes[0].plot([0,1],[es[m/2-1],es[m/2-1]],color='blue',lw=2,zorder=2)
        axes[0].plot([0,1],[es[m/2],es[m/2]],color='green',lw=2,zorder=2)
        axes[0].plot([0,1],np.vstack([es[m/2+1:],es[m/2+1:]]),color='black',lw=1,zorder=1)
        axes[1].plot(range(m),np.abs(evs[:,m/2-1]),label='$E$=%s'%decimaltostr(es[m/2-1]),color='b')
        axes[1].plot(range(m),np.abs(evs[:,m/2]),label='$E$=%s'%decimaltostr(es[m/2]),color='g')
        plt.legend()
        pdb.set_trace()
        plt.savefig('%s/%s_%s.png'%(tba.dout,tba,'EDGE'))
        plt.close()

def edtasks(parameters,lattice,terms,jobs=('GSE',)):
    import HamiltonianPy.ED as ED
    import HamiltonianPy.Beta.TrED as TrED
    ns,ne=len(lattice),len(lattice)/2-(1 if len(lattice)%2==0 and len(lattice.vectors)==0 else 0)
    if 'GSE' in jobs:
        with open('result/ed/%s_%s_%s_GSE_ED.dat'%(name,lattice.name,'_'.join('%s'%decimaltostr(value,n=10) for value in parameters.values())),'w') as fout:
            for i in xrange(ne+1):
                basis=FBasis(ns*2,ne,ne*0.5-i)
                ed=edconstruct(parameters,basis,lattice,terms)
                ed.register(ED.EIGS(name='GSE',ne=1,run=ED.EDEIGS))
                fout.write('%-6s %s\n'%((2*i-ne)/2,ed.records['GSE'][0]))
                fout.flush()
    if 'EDGE' in jobs:
        basis,segments=FBasis(ns*2,ne,ne*0.5-1),np.linspace(1.00,1.32,9)
        ed=edconstruct(parameters,basis,lattice,terms)
        ed.register(ED.EL(name='EL',path=BaseSpace(('sd',segments),('dt',segments**2-2)),ns=ns*2,legend=False,run=ED.EDEL))
        ed.summary()
    if 'MGN' in jobs:
        basis=TrED.TrFBasis(FBasis(ns*2,ne,ne*0.5-1),dk=4,nk=ns/2)
        ed=tredconstruct(parameters,basis,lattice,terms)
        ed.register(TrED.EB(name='TREDEB',ns=2,run=TrED.TrFEDEB))
        ed.summary()

def edgap(fname,parameters,lattice,terms):
    import HamiltonianPy.Beta.TrED as TrED
    assert len(lattice.vectors)==1 and len(lattice)%2==0
    with open('result/ed/%s_GAP_ED.dat'%fname,'w') as fout:
        basis=TrED.TrFBasis(FBasis(len(lattice)*2,len(lattice)/2,0.5*len(lattice)/2-1),dk=4,nk=len(lattice)/2)
        for us in np.linspace(0.02,0.78,39):
            parameters['Us']=us
            ed=tredconstruct(parameters,basis,lattice,terms)
            ed.register(TrED.EB(name='TREDEB',ns=2,savedata=False,plot=False,run=TrED.TrFEDEB))
            result=ed.records['TREDEB']
            fout.write('%1.2f %s %s\n'%(us,result[:,1].min()-result[:,2].max(),result[0,1]-result[0,2]))
            fout.flush()

def fbfmtasks(parameters,lattice,terms,interactions,nk=800,jobs=('EB',)):
    import HamiltonianPy.FBFM as FB
    ns,ne=len(lattice),len(lattice)/2-(1 if len(lattice)%2==0 and len(lattice.vectors)==0 else 0)
    basis=FB.FBFMBasis(BZ=FBZ(lattice.reciprocals,nks=(nk,)) if len(lattice.vectors)>0 else None,filling=Fraction(ne,ns*2))
    fbfm=fbfmconstruct(parameters,basis,lattice,terms,interactions)
    if 'EB' in jobs:
        fbfm.register(FB.EB(name='EB%s'%nk,path='L:G1-G2' if len(lattice.vectors)>0 else None,ne=2 if len(lattice.vectors)>0 else ne*2,run=FB.FBFMEB))
        fbfm.summary()
    if 'BP' in jobs:
        fbfm.register(BP(name='BP',path='L:G1-G2',ns=(0,1,2,3),run=FB.FBFMBP))
        fbfm.summary()
    if 'POS' in jobs:
        fbfm.register(POS(name='POS',ns=[0]+[ns/2+count for count in xrange(-2,2)],run=FB.FBFMPOS))
        fbfm.summary()
    if 'EDGE' in jobs:
        fbfm.register(FB.EB(name='EDGE',path=BaseSpace(('Us',np.linspace(0.02,0.78,39))),ne=ns*2,run=FB.FBFMEB))
        fbfm.summary()

def fbfmVtasks(parameters,nk=800):
    from scipy.linalg import eigh
    import matplotlib.pyplot as plt
    assert parameters['dt']==parameters['sd']**2-2
    engine='%s_%s_EB%s'%(name,'_'.join(decimaltostr(value) for value in parameters.itervalues()),nk)
    result=np.zeros((nk,nk+1))
    print '%s: '%nk,
    for i in xrange(nk):
        print '%s%s'%(i,',' if i<nk-1 else ''),
        result[i,0]=i
        result[i,1:]=eigh(fbfmmatrix(parameters['sd'],parameters['Us'],parameters['Ud'],parameters['V'],i,nk),eigvals_only=True)
    np.savetxt('result/fbfm/%s.dat'%engine,result[:,:nk/2+1])
    plt.plot(result[:,0],result[:,1:nk/2+1])
    plt.title(engine)
    plt.pause(1.0)
    plt.savefig('result/fbfm/%s.png'%engine)
    plt.close()

def fbfmgap(fname,parameters,lattice,terms,interactions):
    import HamiltonianPy.FBFM as FB
    assert len(lattice.vectors)==1 and len(lattice)==2
    with open('result/fbfm/%s_GAP_FBFM_EX_EX.dat'%fname,'w') as fout:
        for us in np.linspace(0.409,0.411,2):
            parameters['Us']=us
            for i,n in enumerate((800,)):
                basis=FB.FBFMBasis(BZ=FBZ(lattice.reciprocals,nks=(n,)),filling=Fraction(1,4))
                fbfm=fbfmconstruct(parameters,basis,lattice,terms,interactions)
                fbfm.log<<'%s\n'%fbfm
                fbfm.register(FB.EB(name='EB',path='L:G1-G2',ne=2,savedata=False,plot=False,run=FB.FBFMEB))
                result=fbfm.records['EB']
                fout.write('%s %s %s'%('%1.3f'%us if i==0 else '',result[:,2].min()-result[:,1].max(),result[n/2,2]-result[n/2,1]))
                fout.flush()
            fout.write('\n')

def fbfmbp(fname,parameters,lattice,terms,interactions):
    import HamiltonianPy.FBFM as FB
    assert len(lattice.vectors)==1 and len(lattice)==2
    with open('result/fbfm/%s_BP_FBFM.dat'%fname,'w') as fout:
        for us in np.linspace(0.02,0.78,39):
            parameters['Us']=us
            basis=FB.FBFMBasis(BZ=FBZ(lattice.reciprocals,nks=(800,)),filling=Fraction(1,4))
            fbfm=fbfmconstruct(parameters,basis,lattice,terms,interactions)
            fbfm.register(BP(name='BP',path='L:G1-G2',ns=(0,1,2,3),savedata=False,plot=False,run=FB.FBFMBP))
            fout.write('%1.2f %s\n'%(us,' '.join(str(bp) for bp in fbfm.apps['BP'].bps)))
            fout.flush()

if __name__=='__main__':
    mkl.set_num_threads(1)
    Engine.DEBUG=True

    m=8
    factor=1.25
    delta=0.0
    parameters=OrderedDict()
    parameters['t']=1.0
    parameters['sd']=factor
    parameters['dt']=factor**2-2+delta

    # tba
    #tbatasks(parameters,S2x('1P-1O',nneighbour),[t,sd,dt],jobs=('EB',))
    #tbatasks(parameters,S1('%sO-1O'%(m+0),nneighbour),[t,sd,dt],tbaidfmap=tbaidfa,jobs=('EDGE',))
    #tbatasks(parameters,S1('%sO-1O'%(m+1),nneighbour),[t,sd,dt],tbaidfmap=tbaidfa,jobs=('EDGE',))
    #tbatasks(parameters,S1('%sO-1O'%(m+1),nneighbour),[t,sd,dt],tbaidfmap=tbaidfb,jobs=('EDGE',))

    # ed
    parameters['Us']=0.9
    parameters['Ud']=1.0

    #edtasks(parameters,S2x('%sP-1O'%m,nneighbour),[t,sd,dt,Us,Ud],jobs=('GSE',))
    #edtasks(parameters,S2x('%sP-1O'%m,nneighbour),[t,sd,dt,Us,Ud],jobs=('MGN',))
    #edtasks(parameters,S2x('%sO-1O'%m,nneighbour),[t,sd,dt,Us,Ud],jobs=('GSE',))
    #edtasks(parameters,S2x('%sO-1O'%m,nneighbour),[t,sd,dt,Us,Ud],jobs=('EDGE',))
    #edgap('S2x(%s)'%m,parameters,S2x('%sP-1O'%m,nneighbour),[t,sd,dt,Us,Ud])

    # fbfm
    #fbfmtasks(parameters,S2x('1P-1O',nneighbour),[t,sd,dt],[Us,Ud],nk=100,jobs=('EB',))
    #fbfmtasks(parameters,S2x('1P-1O',nneighbour),[t,sd,dt],[Us,Ud],nk=800,jobs=('BP',))
    #fbfmtasks(parameters,S2x('60O-1O',nneighbour),[t,sd,dt],[Us,Ud],jobs=('POS',))
    #fbfmtasks(parameters,S2x('60O-1O',nneighbour),[t,sd,dt],[Us,Ud],jobs=('EDGE',))
    #fbfmgap('S2x',parameters,S2x('1P-1O',nneighbour),[t,sd,dt],[Us,Ud])
    #fbfmbp('S2x',parameters,S2x('1P-1O',nneighbour),[t,sd,dt],[Us,Ud])

    # fbfm with V
    parameters['V']=0.97
    #fbfmVtasks(parameters,nk=100)
