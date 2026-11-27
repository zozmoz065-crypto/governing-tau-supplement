import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams['font.family']='DejaVu Sans'
INK='#1b365d'; BL='#257de0'; OR='#f47230'; AQ='#00a888'; GR='#666666'
fig,ax=plt.subplots(figsize=(10,5.6)); ax.set_xlim(0,100); ax.set_ylim(0,56); ax.axis('off')
def box(x,y,w,h,t,ec=INK,fc='white',fs=8.5,dash=False,bold=False):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.3,rounding_size=1.2',ec=ec,fc=fc,lw=1.4,ls='--' if dash else '-')
    ax.add_patch(p); ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=fs,color='#191919',weight='bold' if bold else 'normal',linespacing=1.3)
    return (x,y,w,h)
def arr(p1,p2,c=INK,dash=False,rad=0.0):
    ax.add_patch(FancyArrowPatch(p1,p2,arrowstyle='-|>',mutation_scale=11,color=c,lw=1.3,ls='--' if dash else '-',connectionstyle=f'arc3,rad={rad}'))
box(1,41,16,9,'Signal volume λ\nSignal quality q\nAbsorptive\ncapacity A',fc='#eaf2fc',fs=8)
box(22,41,17,9,'Triage queue\nμ = μ(R, G)\nτ$^{tri}$ convex in ρ = λ$_e$/μ\n(Results 1–2)',fc='#eaf2fc')
box(45,41,15,9,'Decision\nτ$^{dec}$ = f(c, ζ)\n(Result 3)')
box(66,41,15,9,'Enforcement\nτ$^{enf}$ = f(w, m)\n(Result 3)')
box(86,43,13,6.5,'Safety\nconstraints L$_c$',ec=GR,dash=True,fs=8)
box(3,27,13,7,'Security\nresources R',ec=OR,fc='#fef1ea')
box(22,27,13,7,'Governance\nmaturity G',ec=OR,fc='#fef1ea')
box(1,10,36,9,'Conditional complementarity (Result 1b)\nR, G complements in τ$^{tri}$ only if ρ < 1 − 2σ;\nG alone lowers the floor F(G) = E[τ$^{dec}$] + E[τ$^{enf}$]',ec=OR,fc='#fff8f3',fs=7.8)
box(58,24,22,9,'τ$_{tot}$ = τ$_{pre}$ + τ$^{tri}$ + τ$^{dec}$ + τ$^{enf}$\n(two-clock exposure window)',fc=INK,ec=INK,fs=8.5)
ax.texts[-1].set_color('white')
box(84,24,15,9,'Recurrence\nhazard\n(Result 4)',ec=AQ,fc='#e6f7f3')
box(84,8,15,7,'Detection\nintensity D(G)',ec=GR,dash=True,fs=8)
box(45,6,30,8,'Double-loop rule revision (Layer 3)\nraises a(q), ζ; shortens c (P4)',ec=AQ,fc='#e6f7f3',fs=8)
arr((17.6,45.5),(22,45.5)); arr((39.6,45.5),(45,45.5)); arr((60.6,45.5),(66,45.5)); arr((92.5,42.6),(81.6,46.5),c=GR,dash=True)
arr((73.5,40.5),(73.5,33.4)); arr((80.6,28.5),(84,28.5)); arr((95,15.3),(95,23.6),c=GR,dash=True)
arr((11,34.4),(25,40.5),c=OR); arr((28.5,34.4),(30.5,40.5),c=OR); arr((33,34.4),(50,40.5),c=OR)
arr((88,23.5),(75.5,10),c=AQ,rad=-0.2); arr((44.6,11),(30,26.6),c=AQ)
ax.plot([3,8],[2.5,2.5],color=INK,lw=1.3); ax.text(9,2.5,'hypothesized mechanism',va='center',fontsize=7.5)
ax.plot([33,38],[2.5,2.5],color=GR,lw=1.3,ls='--'); ax.text(39,2.5,'boundary condition / bias control',va='center',fontsize=7.5)
ax.plot([68,73],[2.5,2.5],color=AQ,lw=1.3); ax.text(74,2.5,'learning feedback',va='center',fontsize=7.5)
fig.savefig('fig2.png',dpi=300,bbox_inches='tight'); 
