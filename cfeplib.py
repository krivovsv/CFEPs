def comp_Zca(lx, a, dt=1, strict=False, dx=1e-3, zcmin=1e-8, mindx=1e-3):
    import math
    dzc={}
    tmax=len(lx)
    for i in range(0,tmax-dt):
        x=lx[i+dt]
        lastx=lx[i]
        d=abs(x-lastx)
        if a<0 and d<mindx:d=mindx
        if a==0: d=1 
        else: d=float(d)**a
        if dx!=None and not strict: 
            x=math.floor(x/dx)*dx
            lastx=math.floor(lastx/dx)*dx
        if lastx<x:
            dzc[lastx]=dzc.get(lastx,0)+d
            dzc[x]=dzc.get(x,0)-d
        else:  
            dzc[x]=dzc.get(x,0)+d
            dzc[lastx]=dzc.get(lastx,0)-d
    keys=list(dzc.keys())
    keys.sort()
    lx=[]
    ly=[]
    z=0
    for x in keys:
        lx.append(x)
        ly.append(float(z)/dt/2)
        z=z+dzc[x]
        if strict:
            lx.append(x)
            y=float(z)/dt/2
            if y<zcmin: y=0.0
            ly.append(y)
    if not strict:
        ly[0]=ly[1] # ly[0] equals 0, but a non-zero value is more convenieint for further analysis
    return lx,ly 
    
def comp_diffusion_euler(force,D,dt,nsteps,isave=10,x0=0):
    import random,numpy
    dt0=dt/isave
    x=x0
    traj=numpy.zeros((nsteps),'f')
    for istep in range(nsteps):
        for isteps2 in range(isave):
            sigma=(2*D(x)*dt0)**0.5
            dw=random.gauss(0,sigma)
            x=x+force(x)*dt0+dw
        traj[istep]=x
    return traj

def comp_Zh(lx,dx):
    import math
    zh={}
    for x in lx:    
        x=math.floor(x/dx+0.5)*dx
        zh[x]=zh.get(x,0)+1
    for x in zh:zh[x]=float(zh[x])/dx
    lx=list(zh.keys())
    lx.sort()
    ly=[zh[x] for x in lx]
    return lx,ly
    
def comp_ekn_tp(traj,x0,x1,dx=None,dt=1):
    import math
    def process(traj): # process a TP
        n=len(traj)
        if n<2:return
        for i in range(1,n): # from i-dt to i
            j=i-dt
            if j<0:j=0
            key=traj[j],traj[i]
            ekn[key]=ekn.get(key,0)+1
        for i in range(max(n-dt,1),n-1):
            key=traj[i],traj[-1]
            ekn[key]=ekn.get(key,0)+1
        if dt>n-1:
            key=traj[0],traj[-1]
            ekn[key]=ekn.get(key,0)+dt-n+1
    
    ekn={}
    ok=False
    lx=[]
    if dx!=None: 
        x0=math.floor(x0/dx)*dx
        x1=math.floor(x1/dx)*dx
    for x in traj:
        if dx!=None:x=math.floor(x/dx)*dx
        if x<=x0:
            lx.append(x0)
            if ok:process(lx)
            lx=[x0]
            ok=True
            continue
        if x>=x1:  
            lx.append(x1)
            if ok:process(lx)
            lx=[x1]
            ok=True
            continue
        lx.append(x)
    process(lx)
    for ij in ekn:ekn[ij]=float(ekn[ij])/dt   
    return ekn

def comp_Zca_ekn(ekn,a,dx=None,eps=None,zcmin=1e-8,initial0=False):
    import math
    dzc={}
    for x,y in ekn:
        d=abs(y-x)**a*ekn[(x,y)]
        if dx!=None: 
            x=math.floor(x/dx)*dx
            y=math.floor(y/dx)*dx
        if y<x:
            dzc[y]=dzc.get(y,0)+d
            dzc[x]=dzc.get(x,0)-d
        else:  
            dzc[x]=dzc.get(x,0)+d
            dzc[y]=dzc.get(y,0)-d
    keys=list(dzc.keys())
    keys.sort()
    lx=[]
    ly=[]
    z=0
    for x in keys:
        lx.append(x)
        ly.append(float(z)/2)
        z=z+dzc[x]
        if eps!=None:
            lx.append(x+eps)
            y=float(z)/2
            if y<zcmin: y=0
            ly.append(y)
    if not initial0:ly[0]=ly[1] # ly[0] equals 0, but a non-zero value is convenieint for further analysis
    return lx,ly 
    
def to_committor(traj,dx,x0,x1):
    import math
    lx,lzc1=comp_Zca(traj,a=1,dx=dx)
    x2q={}
    q=0
    for x,zc1 in zip(lx,lzc1):
        if x<=x0 or x>=x1:continue
        x2q[x]=q,1./zc1
        q+=1./zc1
    qm=q
    lq=[]
    for x in traj:
        xi=math.floor(x/dx)*dx
        if xi<=x0:lq.append(0)
        elif xi>=x1:lq.append(1)
        else:
            q,dqdx=x2q[xi]
            lq.append((q+(x-xi)*dqdx)/qm)
    return lq
    
def comp_eval(traj,dt=1):
    from math import log
    xx=0
    xx1=0
    for i in range(len(traj)-dt):
        xx=xx+traj[i]**2
        xx1=xx1+traj[i+dt]*traj[i]
    xx1=xx1/xx
    return -log(xx1)/dt


def comp_Z0c1(lx,dx,strict=False):
    import math
    dzc={}
    tmax=len(lx)
    for i in range(tmax):
        x=lx[i]
        lastx=0
        d=abs(x-lastx)
        if dx!=None and not strict: 
            x=math.floor(x/dx)*dx
            lastx=math.floor(lastx/dx)*dx
        if lastx<x:
            dzc[lastx]=dzc.get(lastx,0)+d
            dzc[x]=dzc.get(x,0)-d
        else:  
            dzc[x]=dzc.get(x,0)+d
            dzc[lastx]=dzc.get(lastx,0)-d
    keys=list(dzc.keys())
    keys.sort()
    lx=[]
    ly=[]
    z=0
    for x in keys:
        lx.append(x)
        ly.append(float(z))
        z=z+dzc[x]
    if not strict:ly[0]=ly[1]
    return lx,ly 

def comp_theta(traj,ldt,tinf,dx):
    from math import exp, log
    mu=comp_eval(traj,tinf)
    lx,lz0c1=comp_Z0c1(traj,dx=dx) # for dt>1 z0c1 needs rescaling as z0c1/dt
    ltheta=[]
    for dt in ldt:
        lx,lzc1=comp_Zca(traj,a=1,dx=dx,dt=dt)
        sc=(1-exp(-mu*dt))/dt # the 1/dt factor is to rescale z0c1 -> z0c1/dt
        ly=[-log(zc1/z0c1/sc) for zc1,z0c1 in zip(lzc1,lz0c1)]
        ltheta.append(ly)
    return lx,ltheta


