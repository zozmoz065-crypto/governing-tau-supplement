import math, json
# Structural latency (Result 3, zeta=0): (c/2 + d0) + (x0 + m*w/2)
def S(c,w,d0=2,x0=2,m=1.0): return c/2+d0+x0+m*w/2
rows=[]
for target,lab in [(174,'All KEV (median 174 d)'),(137,'Critical KEV (median 137 d)')]:
  for c,w in [(7,30),(14,30),(30,60)]:
    for pre in [0,30]:
      Q=target-S(c,w)-pre; M=Q/math.log(2)
      r={'target':lab,'c':c,'w':w,'pre':pre,'S':S(c,w),'Q':Q}
      for s in [1,2,5]: r[f'rho_s{s}']=1-s/M
      rows.append(r)
for r in rows: print({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()})
json.dump(rows,open('calib.json','w'),indent=1)
