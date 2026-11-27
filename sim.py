import numpy as np, math, json, heapq
rng=np.random.default_rng(20261127)
out={}
# ---------- DES: single server, non-preemptive 2-class priority; general service ----------
def des(lam, mu, s, a, n=300000, svc='exp', cv2=1.0, prio=True):
    t=0.0; arr=[]
    inter=rng.exponential(1/lam,n); times=np.cumsum(inter)
    crit=rng.random(n)<s
    hi=crit & (rng.random(n)<a) if prio else np.zeros(n,bool)
    if svc=='exp': S=rng.exponential(1/mu,n)
    else:
        sig2=math.log(1+cv2); m=math.log(1/mu)-sig2/2
        S=rng.lognormal(m,math.sqrt(sig2),n)
    q=[]; i=0; free=0.0; soj=np.empty(n)
    while i<n or q:
        # admit all arrivals up to when server free
        if not q:
            free=max(free,times[i])
        while i<n and times[i]<=free:
            heapq.heappush(q,(0 if hi[i] else 1, times[i], i)); i+=1
        if not q: continue
        _,ta,j=heapq.heappop(q)
        start=max(free,ta); free=start+S[j]; soj[j]=free-ta
    w=slice(n//10,None)  # burn-in
    return soj[w][crit[w]].mean(), soj[w].mean()
def eq7(lam,mu,s,a):
    rho=lam/mu; r1=a*s*lam/mu; W0=lam/mu**2
    return 1/mu + a*W0/(1-r1) + (1-a)*W0/((1-r1)*(1-rho))
def pk(lam,mu,cv2): rho=lam/mu; ES=1/mu; ES2=(1+cv2)*ES**2; return ES+lam*ES2/(2*(1-rho))
def erlangc(lam,mu_each,c):
    A=lam/mu_each; rho=A/c
    s_=sum(A**k/math.factorial(k) for k in range(c)); last=A**c/math.factorial(c)/(1-rho)
    Pw=last/(s_+last); return 1/mu_each + Pw/(c*mu_each-lam)
if __name__=='__main__':
    mu=1.0; s=0.15; a=0.9
    rows=[]
    for rho in [0.5,0.7,0.8,0.9,0.95]:
        lam=rho*mu
        simc,_=des(lam,mu,s,a)
        simf,allf=des(lam,mu,s,a,prio=False)
        rows.append(dict(rho=rho, fcfs_an=1/(mu-lam), fcfs_sim=allf, crit_an=eq7(lam,mu,s,a), crit_sim=simc))
    out['validation']=rows
    # robustness: heavy-tailed service (lognormal CV^2=4) via P-K, and M/M/c Erlang C
    rob=[]
    for rho in [0.5,0.7,0.8,0.9,0.95]:
        rob.append(dict(rho=rho, mm1=1/(1-rho), mg1_cv4=pk(rho,1,4.0), mmc3=erlangc(3*rho,1.0,3)*1.0,
                        mg1_sim=des(rho,1,s,a,svc='logn',cv2=4.0,prio=False)[1]))
    out['robust']=rob
    # priority under lognormal service (Result 2 robustness)
    out['prio_logn']=[dict(rho=r, crit=des(r,1,s,a,svc='logn',cv2=4.0)[0], fcfs=des(r,1,s,a,svc='logn',cv2=4.0,prio=False)[1]) for r in [0.8,0.9,0.95]]
    json.dump(out,open('sim.json','w'),indent=1)
    if 0:
        print(k)
        for r in v: print({kk:round(vv,3) for kk,vv in r.items()})
