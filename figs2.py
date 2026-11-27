import numpy as np, math, json, sys
sys.path.insert(0,'.')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import importlib
C1,C2,C3='#257de0','#f47230','#00b894'; INK='#191919'; INK2='#666666'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.edgecolor':INK2,'axes.labelcolor':INK,
 'xtick.color':INK2,'ytick.color':INK2,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,
 'grid.color':'#e6e6e6','grid.linewidth':0.6,'lines.linewidth':2,'legend.frameon':False})
res={}
res['validation']=json.load(open('figs.json'))['validation']
from sim import eq7
def pk(lam,mu,cv2): rho=lam/mu; ES=1/mu; ES2=(1+cv2)*ES**2; return ES+lam*ES2/(2*(1-rho))
def erlangc(lam,mu_each,c):
    A=lam/mu_each; rho=A/c
    s_=sum(A**k/math.factorial(k) for k in range(c)); last=A**c/math.factorial(c)/(1-rho)
    Pw=last/(s_+last); return 1/mu_each + Pw/(c*mu_each-lam)
# ---------- Figure 4: conditional complementarity ----------
def tau_ces(Rv,Gv,lam,sigma,xi=0.5,k=1.0):
    if abs(sigma-1)<1e-9: m=k*Rv**xi*Gv**(1-xi)
    else:
        r=1-1/sigma; m=k*(xi*Rv**r+(1-xi)*Gv**r)**(1/r)
    return 1/(m-lam) if m>lam else np.nan
def cross(Rv,Gv,lam,sigma,h=1e-3):
    f=lambda a,b: tau_ces(a,b,lam,sigma)
    return (f(Rv+h,Gv+h)-f(Rv+h,Gv-h)-f(Rv-h,Gv+h)+f(Rv-h,Gv-h))/(4*h*h)
# grid check of the analytical boundary rho* = 1-2 sigma (R=G=1 so mu=k=1)
S=np.linspace(0.05,1.0,60); P=np.linspace(0.02,0.97,60); Z=np.zeros((len(P),len(S)))
mism=0;tot=0
for j,sg in enumerate(S):
    for i,rr in enumerate(P):
        z=cross(1.0,1.0,rr,sg); Z[i,j]=np.sign(z)
        if abs(rr-(1-2*sg))>0.02:
            tot+=1; mism+= (z<0)!=(rr<1-2*sg)
res['boundary_check']=dict(points=tot,mismatches=int(mism))
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0))
ax[0].contourf(S,P,Z,levels=[-2,0,2],colors=['#d6e7fa','#fde3d4'])
ax[0].plot(S[S<=0.5],1-2*S[S<=0.5],color=INK,lw=1.5)
ax[0].text(0.08,0.15,'Complements\nin latency\n(∂²τ/∂R∂G < 0)',fontsize=7.5,color=INK)
ax[0].text(0.55,0.6,'Substitutes in latency\n(∂²τ/∂R∂G > 0)',fontsize=7.5,color=INK)
ax[0].text(0.27,0.52,'ρ = 1 − 2σ',fontsize=7.5,color=INK,rotation=-50)
ax[0].axvline(1,color=INK2,lw=0.8,ls=':'); ax[0].text(0.83,0.05,'Cobb–\nDouglas',fontsize=7,color=INK2)
ax[0].set_xlabel('Elasticity of substitution σ between R and G'); ax[0].set_ylabel('Utilization ρ')
ax[0].set_title('(a) Sign of the latency cross-partial',fontsize=9,loc='left'); ax[0].grid(False)
# panel b: structural floor (CD), days
lam=8.0; kap=3.0
def floorG(G):
    zeta=0.15*(G-1); c=14/G; return (1-zeta)*(c/2+1.0)+zeta*0.25 + (1.0+0.5*14/2)
Rg=np.linspace(2,25,300); fl={}
for G,col,ls,lab in [(1.5,C2,'--','Low G = 1.5'),(4.5,C1,'-','High G = 4.5')]:
    t=[]
    for Rv in Rg:
        m=kap*Rv**0.5*G**0.5; t.append(1/(m-lam)+floorG(G) if m>lam*1.001 else np.nan)
    ax[1].plot(Rg,t,color=col,ls=ls,label=lab); ax[1].axhline(floorG(G),color=col,lw=0.9,ls=':')
    fl[G]=floorG(G)
ax[1].set_ylim(0,20); ax[1].set_xlabel('Security resources R (analyst FTE-equivalent)'); ax[1].set_ylabel('E[τ] (days)')
ax[1].set_title('(b) Spending cannot breach the floor F(G)',fontsize=9,loc='left'); ax[1].legend(fontsize=7.5)
ax[1].text(20,fl[1.5]+0.4,'floor F(G=1.5)',fontsize=7,color=INK2); ax[1].text(20,fl[4.5]-1.3,'floor F(G=4.5)',fontsize=7,color=INK2)
fig.tight_layout(); fig.savefig('fig4.png',dpi=300); plt.close()
res['floors']=fl
# ---------- Figure 5: interior optimum + detection paradox ----------
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0))
tg=np.linspace(0,30,400); thLx=1.0; phi0=0.5; s0=5.0; opts={}
for Lc,col,ls,lab in [(20,C1,'-',r'IT office systems ($L_c$ = 20)'),(60,C2,'--',r'Finance core ($L_c$ = 60)'),(200,C3,':',r'Clinical / OT ($L_c$ = 200)')]:
    L=thLx*tg+Lc*phi0*np.exp(-tg/s0); ax[0].plot(tg,L,color=col,ls=ls,label=lab)
    to=s0*math.log(Lc*phi0/(thLx*s0)); opts[Lc]=to
    ax[0].plot([to],[thLx*to+Lc*phi0*math.exp(-to/s0)],'o',color=col,ms=6,mfc='white')
ax[0].set_xlabel('Governance latency τ (days)'); ax[0].set_ylabel(r'Expected loss $\Lambda(\tau)$ (units of $\theta L_x$)')
ax[0].set_title(r'(a) $\tau^{opt}$ rises with change-failure loss',fontsize=9,loc='left'); ax[0].legend(fontsize=7.5); ax[0].set_ylim(0,60)
res['tau_opt']=opts
Gg=np.linspace(1,5,200); th=0.2
D=1/(1+np.exp(-1.5*(Gg-2))); ttot=4+16/Gg
ax[1].plot(Gg,th*ttot,color=C1,label=r'True expected events $\theta\tau_{tot}$')
ax[1].plot(Gg,D*th*ttot,color=C2,ls='--',label=r'Observed events $D(G)\theta\tau_{tot}$')
ipk=np.argmax(D*th*ttot); ax[1].plot([Gg[ipk]],[(D*th*ttot)[ipk]],'o',color=C2,ms=6,mfc='white')
ax[1].set_xlabel('Governance maturity G'); ax[1].set_ylabel('Exploitation events per dyad')
ax[1].set_title('(b) Detection paradox (Result 4)',fontsize=9,loc='left'); ax[1].legend(fontsize=7.5); ax[1].set_ylim(0,4.5)
res['paradox_peak_G']=float(Gg[ipk])
fig.tight_layout(); fig.savefig('fig5.png',dpi=300); plt.close()
json.dump(res,open('figs.json','w'),indent=1,default=float)
print(json.dumps(res,indent=1,default=lambda x: round(float(x),3)))
