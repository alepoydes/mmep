import numpy as np
import numba as nb
import matplotlib.pyplot as plt
from matplotlib import animation, rc
import matplotlib.transforms as mtransforms
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from IPython.display import HTML
import scipy.sparse as sparse
import scipy.sparse.linalg
import functools
import warnings
import re  
import os  
#from scipy.sparse import csc_matrix

nptype=np.double
nbtype=nb.f8

def orts(direction):
    direction=np.asarray(direction,dtype=np.float); direction/=np.sqrt(np.sum(direction*direction))
    a=np.cross(direction,[0,0,1] if direction[2]<0.9 else [0,1,0])
    a/=np.sqrt(np.sum(a*a))
    b=np.cross(direction,a)
    return direction,a,b
    
def spin_conus(position, direction, vertices, normals, N=10, height=1.):
    position=np.asarray(position,dtype=float)
    phi=2*np.pi*np.arange(0,N, dtype=np.float)/N
    width=height/2.
    a,b=np.cos(phi)*width,np.sin(phi)*width
    direction,e,g=orts(direction)
    direction*=height
    for k in range(N):
        v0=position+direction
        v1=position+a[k-1]*e+b[k-1]*g-direction
        v2=position+a[k]*e+b[k]*g-direction
        vertices.append([v0, v1, v2])
        n=np.cross(v0-v1,v0-v2); n/=np.sqrt(np.sum(n*n))
        normals.append((n+1)/2)
    return 


# Мотивация для оптимизации на многообразии
# f(x)=dot(x,A(x))/2+dot(x,B)=max, g(x)=dot(x,x)-1=0, A=A^*
# f(x+dx)=f(x)+dot(dx,A(x)+B)+dot(dx,A(dx))/2
# g(x+dx)=g(x)+2*dot(x,dx)+dot(dx,dx)
# l(x,k)=f(x)-k*g(x)
# l(x+dx,k)=l(x,k)+[dot(dx,A(x)+B)-k*2*dot(x,dx)]+[dot(dx,A(dx))/2-k*dot(dx,dx)]
# Выбираем a таким, что линейная по dx часть перпендикулярна к x (нормали к поверхности уровня g).
# d_x l(x,a)=A(x)+B-2*k*x
# dot(x,A(x)+B-2*k*x)=0 => k=dot(x,A(x)+B)/dot(x,x)/2
# Таким образом для g(x)=0 => k=dot(x,A(x)+B)/2 и следовательно
# для фиксированных x, таких что g(x)=0:
# L(x+dx)=l(x+dx,k)=l(x,k)+[dot(dx,A(x)+B)-dot(x,A(x)+B)*dot(x,dx)]
#                  +[dot(dx,A(dx))-dot(x,A(x)+B)*dot(dx,dx)]/2
# При фиксированном x функцию L можно рассматривать как квадратичную форму от dx:
# L(x+dx)=dot(dx,A(x,dx))/2+dot(dx,B(x)),
# A(x,dx)=A(dx)-dot(x,A(x)+B)*dx,
# B(x)=A(x)+B-dot(x,A(x)+B)*x.

# Минимизация бисопряженным градиентом на многообразии:
# A(x,dx) - квадратичная часть функционала в окр. точки x на приращении dx,
# A: M x T_x -> T_x
# Β - линейная часть функционалав в окр. точки x,
# B : M -> T_x
# P(x) - проецирует точку x на ближайшую точку на многообразии,
# P : R^n -> M
# Τ(x,y,dx) - паралельно переносит вектор dx из точки x в точку y.
# T : M x M x T_x -> T_y
# Минимизируемая функция
# f(x+dx)=f(x)+dot(dx,A(x,dx))/2-dot(dx,B(x))
# A линейно по dx, следовательно Α(x,0)=0.
def bigc_manifold(A, B, P, T, x0, maxiter=10, epsilon=1e-14, maxstep=np.inf, debug=0):
    def dot(m1,m2): return np.sum(m1*m2)
    x0=P(x0)
    i0=r0=-B(x0) # невязка (лежит в касательном пространстве) в точке x0
    rho0=alpha=omega0=1 # константа
    v0=p0=0 # сопряженные направления (лежат в касательном пространстве)
    for it in range(maxiter):
        # debug
        r0=-B(x0)
        if debug>0 and it%debug==0:
            print(it,':',np.sqrt(dot(r0,r0)))
        #print(it,':',omega0)
        rho1=dot(i0, r0)
        if np.abs(omega0)<epsilon: 
            #print('omega',omega0)
            return x0
        beta=(rho1/rho0)*(alpha/omega0)
        p0=r0+beta*(p0-omega0*v0)
        v0=A(x0,p0)
        iv=dot(i0,v0)
        if np.abs(iv)<epsilon: 
            #print('iv',iv)
            return x0
        alpha=rho1/iv
        lp=np.sqrt(dot(p0,p0))
        #if alpha*lp>maxstep: alpha=maxstep/lp; print('step reduced')
        h=P(x0+alpha*p0)
        # if h is accurate return h
        s=T(x0,h,r0-alpha*v0)
        i0=T(x0,h,i0); p0=T(x0,h,p0); v0=T(x0,h,v0); 
        t=A(h,s)
        tt=dot(t,t)
        if np.abs(tt)<epsilon: 
            #print('omega0',dot(t,s)/tt)
            return h
        omega0=dot(t,s)/tt
        ls=np.sqrt(dot(s,s))
        if omega0*ls>maxstep: omega0=maxstep/ls; print('step reduced')
        x0=P(h+omega0*s)
        # if x1 is accurate then return x1
        r0=T(h,x0,s-omega0*t)
        # index increment
        rho0=rho1
        i0=T(h,x0,i0); p0=T(h,x0,p0); v0=T(h,x0,v0); 
    return x0

# Operation on vector fields.
def normalize(A):
    """Return unit vectors directed as vectors of A"""
    return A/np.expand_dims(np.sqrt(np.sum(A**2,axis=-1)),-1)

def norm(x): return np.sqrt(np.sum(x**2))

def dot(x,y): return np.sum(x*y)

def pwdot(x,y): return np.sum(x*y,axis=-1)

# x - initial point, y - final point, dx - vector at x.
# Возвращает вектор dx параллельно перенесенный из x в y.
def transport(x,y,dx):
    # x и y должны быть единичной длины
    n=normalize(np.cross(x,y)) # ось вращения
    x2=np.cross(x,n) # вектора n, x2 образуют базис в касательном пр-ве в x
    y2=np.cross(y,n) # вектора n, y2 образуют базис в касательном пр-ве в y
    # раскладываем dx по базису в x
    b=dot(dx,n); c=dot(dx,x2); 
    # паралельно перенесенный вектор имеет в базисе y,n,y2 те же координаты,
    # что и x в базисе x,n,x2
    return b*n+c*y2 

def spherical_triangle_area(a,b,c):
    return 2*np.arctan(pwdot(a,np.cross(b,c))/(1+pwdot(a,b)+pwdot(a,c)+pwdot(b,c)))

def project_to_tangent_space(state, grad): 
    pgrad=grad.copy()
    l=np.sum(state**2, axis=-1)
    mask=l<1e-5
    pgrad-=np.expand_dims(np.sum(state*grad, axis=-1),-1)*state
    pgrad[mask]=0
    return pgrad
 
def spher_distance(a,b,flat=False):
    if flat: return np.sqrt(np.sum((a-b)**2))
    else: 
        phi=np.sum(a*b,axis=-1)
        phi[phi<-1]=-1; phi[phi>1]=1
        return np.sqrt(np.sum(np.arccos(phi)**2))

def path_length(path,flat=False):
    return np.sum([spher_distance(a,b,flat=flat) for a,b in zip(path[1:],path[:-1])])

#def projector(state):
#    state=state.reshape((-1,3))
#    n=state.shape[0]
#    proj=np.zeros((n,3,n,3))
#    for k in range(n):
#      proj[k,:,k,:]=np.eye(3)-state[k,:].reshape((3,1))*state[k,:].reshape((1,3))
#    return proj.reshape((3*n,3*n))

def center(point, direction):
    proj=np.array([[[1.,0.],[0.,1.]]])-direction.reshape((-1,2,1))*direction.reshape((-1,1,2))
    mat=proj.sum(axis=0)
    prod=np.sum(np.sum(proj*point.reshape((-1,1,2)),axis=0),axis=-1)
    return np.linalg.solve(mat.astype(np.float), prod.astype(np.float))

def distance(point, origin):
    return np.sqrt(np.sum((point.reshape((-1,point.shape[-1]))-origin.reshape((-1,point.shape[-1])))**2,axis=-1))

def matrix(oper):
    N, M=oper.shape
    mat=np.empty((N, M))
    for n in range(M):
        S=np.zeros(M)
        S[n]=1
        mat[:,n]=oper(S)
    return mat  

# Auxilliary functions
def mesh_array(mesh):
    """Given list of ranges return array of indices of all vetices."""
    return np.concatenate(tuple(map(lambda i: np.expand_dims(i[:],axis=-1), mesh)), axis=-1)

cdict={'red': ((0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 1.0, 1.0)),
        'green': ((0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)),
        'blue': ((0.0, 1.0, 1.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0))
        }
plt.register_cmap(name='blueblackred', data=cdict)

def animate(init, update, N, interval=50):
    fig,data=init()
    anim=animation.FuncAnimation(fig, lambda i: update(i, data), frames=N, interval=interval, blit=False)
    plt.close(anim._fig)
    return HTML(anim.to_html5_video())

@nb.jit(nbtype[:](nbtype[:],nbtype[:]), nopython=True, cache=True)
def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]])

def energy_contributions_names():   
    names=['Anisotropy','Zeeman','Heisenberg','DMI','Dipolar','Total']
    return names

@nb.jit(nb.void(nb.i4[:], nb.i4, nb.i4), nopython=True, cache=True)
def makeind(idx, delta, bc):
    if bc==1:
        for k in range(idx.shape[0]):
            idx[k]=(k-delta)%idx.shape[0]
    else:
        for k in range(idx.shape[0]):
            k2=k-delta
            idx[k]=k2 if k2>=0 and k2<idx.shape[0] else -1

def distinct_spins(a,b,tol=1e-6):
    return np.sum((a-b)**2,axis=-1)>tol*tol

# Magnetic crystalls
class Lattice(object):
    """
    Contains and parameters of magnetic crystal.
    The magnetic crystal is a crystal such that every atom has a magnetic moment $M_n\in\mathbb R^3$, where $n$ runs over a periodic lattice $\Lambda$. 
    All vectors $M_n$ are assumed to have the same legnth $\mu$, that is $M_n=\mu S_n$
    for appropriate choice of $S_n$, $\|S_n\|=1$.
    Below we consider only finite two-dimanensional lattices with periodicity group $\mathbb Z_{N_1}\times\mathbb Z_{N_2}$, then $\Lambda=\mathbb Z_{N_1}\times\mathbb Z_{N_2}\times\tilde\Lambda$, where $\tilde\Lambda$ is the set of all atoms in a unit cell.
    Below we will sometimes expand indices $M_n=M_{n_1,n_2,k}$, $n\in\Lambda$, $n_\cdot\in\mathbb Z_\cdot$, $k\in\tilde\Lambda$.
    """
    @property
    def card(self): return len(self.cell)

    @property
    def dim(self): return len(self.size)

    @property
    def shape(self): return self.size+(self.card, 3)

    @property
    def number_of_atoms(self): return np.prod(self.size)*self.card

    def random(self):
        """Returns random vector field on the lattice"""
        return normalize(np.random.randn(*(self.size+(self.card,3)), dtype=nptype))

    def uniform(self,vec):
        """Returns uniform vector field on the lattice directed as vector vec"""
        return np.tile(normalize(np.reshape(vec,(1,)*self.dim+(1,3))),self.size+(self.card,1), dtype=nptype)

    def __str__(self):
        return "size:{}\nboundary conditions: {}\ntranslation vectors:\n{}\nSpins in unit cell:\n{}\nBonds:\n{}\nExternal field:\n{}\nAnisottopy:\n{}\n{}\nExchange constant: {}\nDzyaloshinskii-Moriya vectors:\n{}\nMu: {}\nGamma: {}\n".format(self.size, self.bc, self.translations,self.cell,self.bonds,self.H,self.K,self.K0,self.J,self.D,self.mu,self.gamma)

    def axes(self):
        """Given a lattice returns list of ranges of every index."""
        return map(lambda i: np.arange(i), self.size)

    def translated_axes(self,shift):
        return map(lambda i, s: np.roll(np.arange(i),s), self.size, np.array(shift))

    def cells_o(self):
        """Return coordinates of origins of all cells"""
        idx=mesh_array(np.meshgrid(*self.axes(), indexing='ij'))
        x=np.dot(idx,self.translations)
        return x.reshape(self.size+(self.dim,))

    def atoms(self):
        a=np.expand_dims(self.cells_o(),-2)
        b=np.concatenate(tuple(map(lambda c: a+c,self.cell)),axis=-2)
        return b.reshape(self.size+(self.card,self.dim))

    def plot(self):
        A=self.cells_o()
        a=A.reshape((-1,self.dim))
        for bond in self.bonds:
            ia=np.ix_(*self.translated_axes(np.array(bond[0])))
            c=a+self.cell[bond[2]]
            b=A[ia].reshape((-1,self.dim))
            d=b+self.cell[bond[1]]
            plt.plot((c[:,0],d[:,0]),(c[:,1],d[:,1]),'-r')
        for k in range(self.card):
            b=a+self.cell[k]
            plt.plot(b[:,0],b[:,1],'.')
        plt.axis('equal'); plt.axis('off')
        omn=np.min(a,axis=0); omx=np.max(a,axis=0); 
        amn=np.min(self.cell,axis=0); amx=np.max(self.cell,axis=0); 
        plt.xlim([omn[0]-(omx[0]-omn[0])*0.01-amn[0],omx[0]+(omx[0]-omn[0])*0.01+amx[0]])
        plt.ylim([omn[1]-(omx[1]-omn[1])*0.01-amn[1],omx[1]+(omx[1]-omn[1])*0.01+amx[1]])   

    def plot_components(self, ergy, names=None, dist=None, ax=None, idx=None, colors=['r','m','k','b','g','k'], markers=['s','<','s','^','v','o'], markercolors=[None,None,None,None,None,'w'], k=3, fontsize=None, markersize=None):
        if ax is None: fig, ax=plt.subplots(1,1)
        if names is None: names=self.energy_contributions_names()
        if dist is None: dist=np.linspace(0.,1.,ergy.shape[0])
        if idx is None: idx=filter(lambda i: ergy[0,i]!=0, range(ergy.shape[1])) 
        for i in idx:
            y=ergy[:,i]-ergy[-1,i]
            ax.plot(dist, y, marker=markers[i], color=colors[i], markersize=markersize, markerfacecolor=markercolors[i])
            if y[k]>0:
                ax.text(dist[k], y.max(), names[i], fontsize=fontsize, ha='left', va='bottom')
            else: 
                ax.text(dist[k], y.min(), names[i], fontsize=fontsize, ha='left', va='top')
        ax.set_xlabel('Distance along path (rad)', fontsize=fontsize)
        ax.set_ylabel('Energy (meV)', fontsize=fontsize)


    def plot_stencil(self):
        linex=[]; liney=[]; vecx=[]; vecy=[]
        for b, D in zip(self.bonds, self.D):
            start=self.cell[b[2]]
            end=-np.dot(b[0],self.translations)+self.cell[b[1]]
            linex.append([start[0],end[0]]); liney.append([start[1],end[1]])
            vecx.append(D[0]); vecy.append(D[1]); 
        linex=np.array(linex); liney=np.array(liney); 
        vecx=np.array(vecx); vecy=np.array(vecy); 
        origin=lambda x: 0.75*x[...,1]+0.25*x[...,0]
        plt.axis('equal')
        plt.quiver(origin(linex),origin(liney),vecx,vecy,angles='xy',units='xy',scale=1.)
        plt.plot(linex.T,liney.T)   

    def plot_angles(self, path, rscale=1., ax=None, style='-', lw=1, origin=None, highlight=None, highlightlw=4, threshold=0.5):
        points=self.atoms()
        if origin is None:
            mask=np.logical_and(path[0,:,:,:,2]>-threshold,path[0,:,:,:,2]<threshold)
            origin=center(points[mask], path[0,mask,:2])
        else:
            origin=np.array(origin)
        dist=distance(points, origin)
        distidx=np.argsort(dist)
        dist=dist[distidx]*rscale
        states=[path[k,:,:,:,2].flatten()[distidx] for k in range(path.shape[0])]
        if ax is None: fig, ax=plt.subplots(1,1)
        for state in states:
            ax.plot(dist, state, style, lw=lw)
        if highlight is not None:
            self.plot_angles(highlight, rscale=rscale, ax=ax, style='-k', lw=highlightlw)
        ax.set_ylim([-1,1])
        return dist, states          

    def plot_field(self,S,scale=1,axis=None,fig=None):
        if axis is None:
            fig,axis=plt.subplots()
        axis.set_aspect('equal')
        axis.set_axis_off()
        s=S.reshape((-1,3))
        b=self.atoms().reshape((-1,self.dim))
        plt.set_cmap('blueblackred')
        Q=axis.quiver(b[:,0],b[:,1],s[:,0],s[:,1],s[:,2],angles='xy',scale_units='xy',pivot="mid",scale=scale,clim=(-1,1))
        omn=np.min(b,axis=0); omx=np.max(b,axis=0)
        axis.set_xlim([omn[0]-1,omx[0]+1])
        axis.set_ylim([omn[1]-1,omx[1]+1])
        return (fig,axis,Q)

    def plot_path(self, path, interval=200):
        def init(): 
            fig,ax,Q=self.plot_field(frame(0))
            return (fig, Q)
        def update(i, Q):
            f=path[i].reshape((-1,3))
            Q.set_UVC(f[:,0],f[:,1],f[:,2])
            return Q
        def frame(i): return path[i]
        return animate(init, update, path.shape[0], interval)
  
    def conus_plot(self, field, fig=None, scale=0.5, axis=None):
        if fig is None: fig=plt.figure(figsize=(10,10))
        if axis is None: axis = fig.add_subplot(111, projection='3d')
        s=field.reshape((-1,3))
        b0=self.atoms().reshape((-1,self.dim))
        b=np.zeros((b0.shape[0],3))
        b[:,:self.dim]=b0
        vertices = []
        normals = []
        height=np.min(np.sum(np.asarray(self.translations)**2,axis=-1))*scale
        for j in range(s.shape[0]):
            spin_conus(b[j],s[j],vertices,normals,height=height)
        collection = Poly3DCollection(vertices, linewidths=0)
        collection.set_facecolors(normals)
        axis.add_collection3d(collection)
        border=2*height
        mn=np.min(b,axis=0); mx=np.max(b,axis=0); center=(mx+mn)/2.; size=np.max(mx-mn)
        axis.set_xlim(center[0]-size/2-border,center[0]+size/2+border)
        axis.set_ylim(center[1]-size/2-border,center[1]+size/2+border)
        axis.set_zlim(center[2]-size/2-border,center[2]+size/2+border)
        return fig, axis

    def coffee_plot(self, state, axis=None, fig=None, spin=0, smooth=False, scale=1):
        if axis is None:
            fig,axis=plt.subplots()
        state=state.astype(np.float)
        matrix=np.eye(3)
        matrix[:2,:2]=np.array(self.translations).T*scale
        szx=self.size[0]*matrix[0,0]+self.size[1]*matrix[0,1]
        szy=self.size[0]*matrix[1,0]+self.size[1]*matrix[1,1]
        interpolation='lanczos' if smooth else 'none'
        im=axis.imshow(state[:,:,spin,2].T, clim=(-1,1), interpolation=interpolation, cmap='YlOrBr', origin='lower',clip_on=True)
        x1, x2, y1, y2 = im.get_extent()
        im._image_skew_coordinate = (x2, y1)
        trans_data = mtransforms.Affine2D(matrix) + axis.transData
        im.set_transform(trans_data)
        axis.set_xlim(-matrix[0,0],szx)
        axis.set_ylim(-matrix[1,1],szy)
        axis.set_axis_off()    
        return (fig, axis, im)

    def lambda_hessian(self, radius=None):
        if(self.dim!=2): raise Exception("Only 2D lattices are supported")
        K2=-2.*self.K; mu2=2.*self.mu
        K0=self.K0
        sx=self.size[0]; sy=self.size[1]; card=self.card
        shift, idxd, idxs = zip(*self.bonds); 
        shift=tuple(map(lambda x: tuple(x), shift))
        selfD=tuple(map(lambda x: tuple(x), self.D))
        selfJ=tuple(self.J)
        pos=self.atoms()
        if radius==None:
            radx=sx-1
            rady=sy-1
        else:
            radx=rady=radius
        bc=np.array(self.bc, dtype=np.int32)
        assert(np.all(bc<=1))
        assert(np.all(bc>=0))

        @nb.jit(nbtype[:,:,:,:](nbtype[:,:,:,:]), nopython=True, cache=True)
        def fun(S):
            if((sx,sy,card,3)!=S.shape): raise Exception("Wrong vector field dimensions")
            ix=np.empty(sx, dtype=np.int32)
            iy=np.empty(sy, dtype=np.int32)          
            dE=np.empty((sx,sy,card,3))
            for x in range(sx): 
                for y in range(sy): 
                    for c in range(card):
                        R=dE[x,y,c]; 
                        R[:]=0
                        for n in range(K0.shape[0]):
                            prod=K0[n,0]*S[x,y,c,0]+K0[n,1]*S[x,y,c,1]+K0[n,2]*S[x,y,c,2];
                            R[0]+=K2[n]*K0[n,0]*prod
                            R[1]+=K2[n]*K0[n,1]*prod
                            R[2]+=K2[n]*K0[n,2]*prod
            for sh, d, s, D, J in zip(shift, idxd, idxs, selfD, selfJ):
                dx, dy = sh
                makeind(ix, dx, bc[0])
                makeind(iy, dy, bc[1])
                for x in range(sx): 
                    x2=ix[x]
                    if x2<0: continue                    
                    for y in range(sy): 
                        y2=iy[y]
                        if y2<0: continue
                        Ss=S[x,y,s]; Sd=S[x2,y2,d]
                        Rs=dE[x,y,s]; Rd=dE[x2,y2,d]
                        Rs-=J*Sd
                        Rd-=J*Ss
                        Rs[0]-=Sd[1]*D[2]-Sd[2]*D[1]
                        Rs[1]-=Sd[2]*D[0]-Sd[0]*D[2]
                        Rs[2]-=Sd[0]*D[1]-Sd[1]*D[0]
                        Rd[0]+=Ss[1]*D[2]-Ss[2]*D[1]
                        Rd[1]+=Ss[2]*D[0]-Ss[0]*D[2]
                        Rd[2]+=Ss[0]*D[1]-Ss[1]*D[0]
            if mu2==0.: return dE

            for d in range(card):
               for s in range(card):
                    for dx in range(-radx,radx+1): 
                        if dx>0: ax=0; bx=dx; 
                        else: ax=-dx; bx=0 
                        for dy in range(-rady,rady+1):
                            if dx==0 and dy==0 and d==s: continue
                            if dy>0: ay=0; by=dy; 
                            else: ay=-dy; by=0
                            u=pos[ax,ay,s]-pos[bx,by,d]
                            r=np.sqrt(u[0]**2+u[1]**2)
                            u/=r

                            makeind(ix, dx, bc[0])
                            makeind(iy, dy, bc[1])
                            for x in range(sx):
                                x2=ix[x]
                                if x2<0: continue
                                for y in range(sy):
                                    y2=iy[y]
                                    if y2<0: continue
                                    m=3*(u[0]*S[x,y,s,0]+u[1]*S[x,y,s,1])
                                    dc0=S[x,y,s,0]-u[0]*m
                                    dc1=S[x,y,s,1]-u[1]*m
                                    dc2=S[x,y,s,2]
                                    m=mu2/(r*r*r)
                                    dE[x2,y2,d,0]+=m*dc0
                                    dE[x2,y2,d,1]+=m*dc1
                                    dE[x2,y2,d,2]+=m*dc2
            return dE
        return fun

    def lambda_energy_contributions(self, radius=None):
        if(self.dim!=2): raise Exception("Only 2D lattices are supported")
        K=self.K; mu=self.mu; H=self.H
        K0=self.K0
        sx=self.size[0]; sy=self.size[1]; card=self.card
        shift, idxd, idxs = zip(*self.bonds); 
        shift=tuple(map(lambda x: tuple(x), shift))
        selfD=tuple(map(lambda x: tuple(x), self.D))
        selfJ=tuple(self.J)
        pos=self.atoms()
        if radius==None:
            radx=sx-1
            rady=sy-1
        else:
            radx=rady=radius
        bc=np.array(self.bc, dtype=np.int32)
        assert(np.all(bc<=1))
        assert(np.all(bc>=0))
       
        @nb.jit(nbtype[:](nbtype[:,:,:,:]), nopython=True, cache=True)
        def fun(S):
            if((sx,sy,card,3)!=S.shape): raise Exception("Wrong vector field dimensions")
            ix=np.empty(sx, dtype=np.int32)
            iy=np.empty(sy, dtype=np.int32)          
            HeisenbergExchange=0.
            DzyaloshinkiiMoriyaInteraction=0.
            DipolarCoupling=0.
            ani=np.zeros(K.shape[0])
            zee=np.zeros(3)
            for x in range(sx): 
                for y in range(sy): 
                    for c in range(card):
                        for n in range(K0.shape[0]):
                            ani[n]+=(S[x,y,c,0]*K0[n,0]+S[x,y,c,1]*K0[n,1]+S[x,y,c,2]*K0[n,2])**2
                        zee+=S[x,y,c]
            Anisotropy=-np.sum(ani*K)
            Zeeman=-(H[0]*zee[0]+H[1]*zee[1]+H[2]*zee[2])

            for sh, d, s, D, J in zip(shift, idxd, idxs, selfD, selfJ):
                dx, dy = sh
                makeind(ix, dx, bc[0])
                makeind(iy, dy, bc[1])
                dmix=0.; dmiy=0.; dmiz=0.; he=0.
                for x in range(sx): 
                    x2=ix[x]
                    if x2<0: continue                    
                    for y in range(sy): 
                        y2=iy[y]
                        if y2<0: continue
                        Ss=S[x,y,s]; Sd=S[x2,y2,d]
                        dmix+=Ss[1]*Sd[2]-Ss[2]*Sd[1]
                        dmiy+=Ss[2]*Sd[0]-Ss[0]*Sd[2]
                        dmiz+=Ss[0]*Sd[1]-Ss[1]*Sd[0]
                        he+=Ss[0]*Sd[0]+Ss[1]*Sd[1]+Ss[2]*Sd[2]
                HeisenbergExchange-=J*he 
                DzyaloshinkiiMoriyaInteraction-=D[0]*dmix+D[1]*dmiy+D[2]*dmiz

            if mu!=0.: 
                for d in range(card):
                    for s in range(card):
                        for dx in range(-radx,radx+1): 
                            if dx>0: ax=0; bx=dx; 
                            else: ax=-dx; bx=0 
                            for dy in range(-rady,rady+1):
                                if dx==0 and dy==0 and d==s: continue
                                if dy>0: ay=0; by=dy; 
                                else: ay=-dy; by=0
                                u0=pos[ax,ay,s,0]-pos[bx,by,d,0]
                                u1=pos[ax,ay,s,1]-pos[bx,by,d,1]
                                r=np.sqrt(u0**2+u1**2)
                                if r==0.: continue                          
                                u0/=r; u1/=r
                                dc=0.
                                makeind(ix, dx, bc[0])
                                makeind(iy, dy, bc[1])
                                for x in range(sx):
                                    x2=ix[x]
                                    if x2<0: continue
                                    for y in range(sy):
                                        y2=iy[y]
                                        if y2<0: continue
                                        Ss=S[x,y,s]; Sd=S[x2,y2,d] 
                                        dc+=3*(Ss[0]*u0+Ss[1]*u1)*(Sd[0]*u0+Sd[1]*u1)-(Sd[0]*Ss[0]+Sd[1]*Ss[1]+Sd[2]*Ss[2])
                                DipolarCoupling+=dc/(r*r*r)
                DipolarCoupling*=-mu
            TotalEnergy=Anisotropy+HeisenbergExchange+DzyaloshinkiiMoriyaInteraction+DipolarCoupling+Zeeman
            result=np.array([Anisotropy,Zeeman,HeisenbergExchange,DzyaloshinkiiMoriyaInteraction,DipolarCoupling,TotalEnergy])
            return result
        return fun

    def lambda_energy(self):
        if(self.H[0]!=0. or self.H[1]!=0.): 
            raise Exception("Magnetic field must be parallel to Oz axis")
        sx=self.size[0]; sy=self.size[1]; card=self.card
        H=self.H[2]
        @nb.jit((nbtype[:,:,:,:],nbtype[:,:,:,:]), nopython=True, cache=True)
        def fun(S, hess):
            E2=0.; E1=0.
            for x in range(sx):
                for y in range(sy):
                    for c in range(card):
                        s=S[x,y,c]; h=hess[x,y,c]
                        E2+=s[0]*h[0]+s[1]*h[1]+s[2]*h[2]
                        E1-=s[2]
                        #E-=(s[0]*h[0]+s[1]*h[1]+s[2]*h[2])/2.
                        hess[x,y,c,2]-=H
                        #E+=(s[0]*h[0]+s[1]*h[1]+s[2]*h[2])
            return (E1*H+E2/2., hess)
        return fun

    def energy(self,S,hess):
        """Compute energy E and energy gradient for the <system> on the state S."""
        dE=hess.reshape((-1,3))
        s=S.reshape((-1,3))
        E=-np.sum(dE*s)/2
        dE-=self.H.reshape((1,3))
        E+=np.sum(dE*s)
        return (E,dE.reshape(S.shape))
    
    def gradient(self,S,hess):
        """Compute energy gradient for the <system> on the state S."""
        return hess-self.H.reshape((1,)*self.dim+(1,3))

    def lambda_restricted_hessian(self, S, hessian=None, mask=None, grad=None, energy=None):
        shp=S.shape[:-1]; shp0=S.shape
        if hessian is None:
            #print('Compiling Hessian')
            hessian=self.lambda_hessian()
        #print('Mark active spins')
        l=np.sqrt(np.sum(S**2, axis=-1)) 
        if not mask is None: mask=np.logical_and(mask, l!=0)
        else: mask=l!=0
        #S/=np.expand_dims(l,-1)
        #S=S.copy(); S[~mask]=0
        idx=np.nonzero(mask)
        N=idx[0].shape[0]
        #print('compute Hessian on the state')
        hess=hessian(S)
        #print('compute gradient and Lagrange multipliers')
        if grad is None and energy is None: 
            energy, grad=self.energy(S, hess)
        elif energy is None: 
            energy=self.energy(S, hess)[0]
        elif grad is None:
            grad=self.gradient(S, hess)
        MU=np.expand_dims(np.sum(grad*S,axis=-1), -1)
        #print('compute basis in the tangent space')
        e0=np.empty(shp0); e0[...,0]=-S[...,1]; e0[...,1]=S[...,0]; e0[...,2]=S[...,2]
        e1=np.cross(e0, S); 
        l=np.sum(e1**2,axis=-1)==0
        e0[l,0]=S[l,0]; e0[l,1]=-S[l,2]; e0[l,2]=S[l,1]
        e1=np.cross(e0, S)
        l=np.sqrt(np.sum(e1**2,axis=-1))
        #assert(N==np.sum(l!=0))
        l[l==0]=1.
        e1/=np.expand_dims(l,-1)
        e2=np.cross(e1, S)
        e1=e1[idx]; e2=e2[idx]

        def embed(P):
            P=P.reshape((-1,2))
            S=np.zeros(shp0)
            S[idx]=P[:,0:1]*e1+P[:,1:]*e2
            return S

        def project(S2):
            S=S2[idx]
            Q=np.empty((N,2))
            Q[:,0]=np.sum(S*e1, axis=-1)
            Q[:,1]=np.sum(S*e2, axis=-1)
            return Q.flatten()

        def restricted_hessian(P):
            #print('embed tangent to N*3d space')
            S=embed(P)          
            #print('compute Hessian of Lagrangian')
            S2=hessian(S)-MU*S
            #print('Project result to the tangent')
            return project(S2)

        oper=sparse.linalg.LinearOperator((2*N, 2*N), restricted_hessian)
        return oper, 2*N, embed, project, energy, grad


    def fourier(self, S):
        return np.fft.fftn(S, axes=np.arange(self.dim), norm="ortho")

    def ifourier(self, S):
        return np.real(np.fft.ifftn(S, axes=np.arange(self.dim), norm="ortho"))

    def translate_fourier(self,vec,S):
        x=map(lambda v, s: np.exp(-2j*np.pi*v*np.fft.fftfreq(s)), vec, self.size)
        x=np.ix_(*x)
        m=functools.reduce(lambda acc, x: acc*x, x, 1)
        return S*m.reshape(self.size+(1,1))

    def translate(self,vec,S):
        return self.ifourier(self.translate_fourier(vec,self.fourier(S)))

    def generator_fourier(self,vec,S):
        x=map(lambda v, s: -2j*np.pi*v*np.fft.fftfreq(s), vec, self.size)
        x=np.ix_(*x)
        m=functools.reduce(lambda acc, x: acc+x, x, 0)
        return S*m.reshape(self.size+(1,1))

    def igenerator_fourier(self,vec,S):
        x=map(lambda v, s: -2j*np.pi*v*np.fft.fftfreq(s), vec, self.size)
        x=np.ix_(*x)
        m=functools.reduce(lambda acc, x: acc+x, x, 0)
        R=S/m.reshape(self.size+(1,1))
        R[np.isnan(R)]=0;
        return R

    def generator(self,vec,S):
        return self.ifourier(self.generator_fourier(vec,self.fourier(S)))

    def restricted_harmonic(self, S, hessian=None, threshold=1e-3, mask=None, grad=None, energy=None):
        oper, N, embed, project, ergy, grad=self.lambda_restricted_hessian(S, hessian=hessian, mask=mask, grad=grad, energy=energy)
        hes=matrix(oper)
        ei, ev=np.linalg.eigh(hes)
        idx=np.argsort(ei)
        ei=ei[idx]
        ev=ev[:,idx]
        # extract zero modes
        msk=np.abs(ei)<threshold
        if np.any(msk):
            zei=ei[msk]
            zev=ev[:,msk]
            msk=np.logical_not(msk)
            ei=ei[msk]
            ev=ev[:,msk]
        else:
            zei=None
            zev=None
        # check if minimum
        msk=ei<0
        if np.any(msk):
            nei=ei[msk]; nev=ev[:,msk]
            msk=np.logical_not(msk)
            pei=ei[msk]; pev=ev[:,msk]
            c=np.array([np.dot(project(np.cross(embed(nev[:,n]),S)), pev) for n in range(nev.shape[1])])
        else:
            pei=ei
            nei=None
            c=None
        return (ergy, pei, nei, c, zei, zev)

    def tune_saddle(self, state, hessian=None, mask=None, maxiter=500, epsilon=1e-8, maxstep=1e-1, degree=2):
        self.tune(state, nnev=1, hessian=hessian, mask=mask, maxiter=maxiter, epsilon=epsilon, maxstep=maxstep, degree=degree)

    def tune(self, state, nnev=0, hessian=None, mask=None, maxiter=500, epsilon=1e-8, maxstep=1e-1, degree=2, inv_mult=0.1, debug=0, zerotol=1e-4):
        if hessian is None:
            hessian=self.lambda_hessian()
        if inv_mult<0: inv_mult=-inv_mult
        ergy0=None
        for it in range(maxiter):
            oper, N, embed, project, ergy, grad=self.lambda_restricted_hessian(state, hessian=hessian, mask=mask)
            if ergy0 is None: ergy0=ergy
            vals, vecs=sparse.linalg.eigsh(oper, k=degree, which='SA')
            grad=project_to_tangent_space(state, grad)
            res=norm(grad)
            if res<epsilon and np.sum(vals<-zerotol)==nnev: break
            projected_grad=project(grad)
            grad_in_subbasis=np.dot(projected_grad, vecs)
            projected_grad-=np.dot(vecs, grad_in_subbasis)
            #print(': {:.3} {}'.format(miser.norm(projected_grad), grad_in_subbasis))
            msk=vals<-zerotol; msk[:nnev]=np.logical_not(msk[:nnev])
            newvals=vals.copy(); newvals[msk]*=-inv_mult
            newvals[np.abs(vals)<zerotol]=np.Infinity;
            step=dot(projected_grad,projected_grad)/dot(projected_grad,oper(projected_grad))
            delta=embed(np.dot(vecs, grad_in_subbasis/newvals)+step*projected_grad)
            l=norm(delta)
            if l>maxstep: delta*=maxstep/l
            state=normalize(state-delta)
            if debug>0 and it % debug==0:
                print('{}: E={:.2} grad={:.2} eigs={} delta={:.2}'.format(it, ergy-ergy0, res, vals,l))
        return state

    def tune_bicg(self, state, hessian=None, maxiter=500, epsilon=1e-8, maxstep=1e-1, debug=0):
        if hessian is None: hessian=self.lambda_hessian()
        def a(x,dx): 
            adx=hessian(dx)
            grad=self.gradient(x,hessian(x))
            l=np.sum(grad*x,axis=-1)
            return adx-l[...,None]*dx
        def b(x):
            grad=self.gradient(x,hessian(x))
            l=np.sum(grad*x,axis=-1)
            return grad-l[...,None]*x
        return bigc_manifold(a,b,normalize,transport,state,maxiter=maxiter,maxstep=maxstep, debug=debug)
        

    def top_charge_triangle(self, state, a, b, c):
        def roll(state, shift):
            return np.roll(np.roll(state,shift[1],axis=1),shift[0],axis=0)
        va=roll(state[...,a[1],:],a[0])
        vb=roll(state[...,b[1],:],b[0])
        vc=roll(state[...,c[1],:],c[0])
        return np.sum(spherical_triangle_area(va,vb,vc))/4/np.pi

    def topological_charge(self, state, k=0):
        t1=self.top_charge_triangle(state,((0,0),k),((1,0),k),((0,1),k))+self.top_charge_triangle(state,((1,1),k),((0,1),k),((1,0),k))
        t2=self.top_charge_triangle(state,((0,0),k),((1,1),k),((1,0),k))+self.top_charge_triangle(state,((0,0),k),((0,1),k),((1,1),k))           
        return (t1,t2)

def rate(initial, transition, kT=1, gammaovermu=1):
    ergy0, pei0, nei0, a0, zei0, zev0=initial
    ergy1, pei1, nei1, a1, zei1, zev1=transition
    if not nei0 is None:
        warnings.warn("Initial state is not a minimum")
    if nei1 is None:    
        warnings.warn("Transition state is not a saddle point")
    if nei1.shape[0]!=1:
        warnings.warn("Transition state must be of first order")
    if ergy0>ergy1:
        warnings.warn("Energy of initial state larger than of transition state")
    exp=np.exp(-(ergy1-ergy0)/kT)
    n=min(pei0.shape[0],pei1.shape[0])
    detr=np.sqrt(np.prod(pei0[:n]/pei1[:n])*np.prod(pei0[n:])/np.prod(pei1[n:]))
    qin=np.sqrt(np.sum(a1*a1*pei1))
    nzei0=0 if zei0 is None else zei0.shape[0]
    nzei1=0 if zei1 is None else zei1.shape[0]
    zeromodes=(2*np.pi*kT)**((nzei0-nzei1)/2.)
    return np.array([gammaovermu/2/np.pi, zeromodes, detr, qin, exp])

class LatticeRohart(Lattice):
    """Geometry and parameters of magnetic crystal with triagonal lattice"""
    def __init__(self, size=(30,30), H=0., K=0., J=1., D=0., gamma=1, mu=0.,bc=(0,0)):
        if(len(size)!=2): raise Exception("Lattice should be 2D")
        # number of cells along each axis
        self.size=size
        # coordinates of atoms in the unit cell
        self.cell=[np.array([0.,0]), np.array([np.cos(np.pi/3),np.sin(np.pi/3)])]
        # coordinates of translation vectors 
        self.translations=[np.array([1,0]),np.array([0,2*np.sin(np.pi/3)])]
        # list of nearest neighbours in the format (<cell>,atom in the cell,atom in the reference cell)
        self.bonds=[([-1,0],0,0),([-1,0],1,1),([0,0],0,1),([0,-1],0,1),([-1,0],0,1),([-1,-1],0,1),]
        # external magentic field
        self.H=np.array([0,0,H])
        # anisotropy value
        self.K=np.array([K])
        # unit anisotropy vector
        self.K0=np.array([[0,0,1]])
        # exchange interaction coefficient for every pair from <bonds>
        self.J=[J for b in self.bonds]
        # D.-M. vectors for every pair from <bonds>
        self.D=[np.cross(normalize(-np.dot(b[0],self.translations)+self.cell[b[1]]-self.cell[b[2]]),[0,0,D]) for b in self.bonds]
        # mu_0*\|S\|^2/8/pi
        self.mu=mu
        # gyromagetic ratio
        self.gamma=gamma
        self.bc=bc  

class LatticeTriagonal(Lattice):
    """Geometry and parameters of magnetic crystal with triagonal lattice"""
    def __init__(self, size=(30,30), H=0., K=0., J=1., D=0., gamma=1., mu=0., bc=(0,0)):
        if(len(size)!=2): raise Exception("Lattice should be 2D")
        # number of cells along each axis
        self.size=size
        # coordinates of atoms in the unit cell
        self.cell=[np.array([0.,0])]
        # coordinates of translation vectors 
        self.translations=[np.array([1,0]),np.array([np.cos(np.pi/3),np.sin(np.pi/3)])]
        # list of nearest neighbours in the format (<cell>,atom in the cell,atom in the reference cell)
        self.bonds=[([1,0],0,0),([0,1],0,0),([-1,1],0,0)]
        # external magentic field
        self.H=np.array([0,0,H])
        # anisotropy value
        self.K=np.array([K])
        # unit anisotropy vector
        self.K0=np.array([[0,0,1]])
        # exchange interaction coefficient for every pair from <bonds>
        self.J=[J,J,J]
        # D.-M. vectors for every pair from <bonds>
        self.D=[D*np.array([0,-1,0]),D*np.array([np.cos(np.pi/6),-np.sin(np.pi/6),0]),D*np.array([np.cos(np.pi/6),np.sin(np.pi/6),0])]
        # gyromagetic ratio
        self.gamma=gamma
                # mu_0*\|S\|^2/8/pi
        self.mu=mu
        self.bc=bc

class LatticeQuadratic(Lattice):
    """Geometry and parameters of magnetic crystal with quadratic lattice"""
    def __init__(self, size=(30,30), H=0., K=0., J=1., D=0., gamma=1., mu=0., bc=(0,0)):
        if(len(size)!=2): raise Exception("Lattice should be 2D")
        self.size=size
        self.cell=[np.array([0.,0])]
        self.translations=np.array([[1,0],[0.,1]])
        self.bonds=[([1,0],0,0),([0,1],0,0)]
        self.H=np.array([0,0,H])
        self.K=np.array([K])
        self.K0=np.array([[0,0,1]])
        self.J=[J,J]
        self.D=[D*np.array([0,-1,0]),D*np.array([1,0,0])]
        self.gamma=gamma  
        self.mu=mu
        self.bc=bc        

class LatticeHexagonal(Lattice):
    """Geometry and parameters of magnetic crystal with hexagonal lattice"""
    def __init__(self, size=(30,30), H=0., K=0., J=1., D=0., gamma=1., mu=0.,bc=(0,0)):
        if(len(size)!=2): raise Exception("Lattice should be 2D")
        self.size=size
        self.cell=[np.array([0.,0]),np.array([0, 1])]
        self.translations=np.array([[1.7320508, 0],[0.8660254, 1.5]])
        self.bonds=[([0,0],1,0),([0,1],1,0),([-1,1],1,0)]
        self.H=np.array([0,0,H])
        self.K=np.array([K])
        self.K0=np.array([[0,0,1]])
        self.J=[J,J,J]
        self.D=D*np.array([[0,1,0],[0.8660254,-0.5,0],[-0.8660254,-0.5,0]])
        self.gamma=gamma           
        self.mu=mu
        self.bc=bc

# Loading data from octave files
class LatticeFromOct(Lattice):
    """Geometry and parameters of magnetic crystal"""
    def __init__(self, octdata, gamma=1):
        self.size=tuple(octdata['SZ'].astype(int).reshape(-1))
        # D - number of symmetry group generators
        D=2 if len(self.size)==2 or self.size[2] else 3 
        self.size=self.size[:D]
        #assert((octdata['BC'].reshape(-1)[:D].astype(int)==1).all())
        self.cell=octdata['CELL'][:,:D]
        self.translations=octdata['TRANSLATIONS'][:D,:D]
        self.bonds=octdata['BONDS'].astype(int)
        self.H=octdata['H'].reshape(-1)
        self.K0=octdata['K0']; #.reshape((-1,3))
        self.K=octdata['K'].reshape(-1)
        self.J=octdata['J'].reshape(-1)
        self.D=octdata['D']
        self.mu=octdata['mu'] if 'mu' in octdata else 0.        
        self.cell=list(map(lambda i: self.cell[i,:D], range(self.cell.shape[0])))
        self.bonds=list(map(lambda i: (self.bonds[i,:D], self.bonds[i,3], self.bonds[i,4]), range(self.bonds.shape[0])))
        self.gamma=gamma 
        self.bc=np.array(octdata['BC'].flatten()[:D],dtype=np.int)

def read_doubles(fileobj, N):
    """Read N lines from fileobj each containing single integer"""
    res=np.empty(N, dtype=nptype)
    for l in range(N): 
        res[l]=np.fromstring(fileobj.readline(), dtype=nptype, sep='\n')
    return res

#def read_sparse(fileobj, nnz, rows, columns):
#    col=np.empty(nnz, dtype=np.int)
#    row=np.empty(nnz, dtype=np.int)
#    data=np.empty(nnz, dtype=np.float64)
#    for l in range(nnz): 
#        dt=map(float, fileobj.readline().split())
#        row[l]=dt[0]
#        col[l]=dt[1]
#        data[l]=dt[2]
#    return csc_matrix((data, (row-1, col-1)), shape=(rows, columns))

def read_array(fileobj):
    s=fileobj.readline()
    if not s: return None
    m=re.match(r"# name: (\w+)", s); assert(m); name=m.group(1)
    #print(name,end=" ")
    m=re.match(r"# type: (.+)", fileobj.readline()); assert(m); tp=m.group(1)
    if tp=='matrix':
        m=re.match(r"# ndims: ([0-9]+)", fileobj.readline()); assert(m); dims=int(m.group(1))
        sz=np.array(list(map(int, fileobj.readline().split())),dtype=np.int)
        assert(len(sz)==dims)
        data=read_doubles(fileobj, sz.prod()).reshape(tuple(sz), order='F')
    elif (tp=='scalar'):
        data=float(fileobj.readline())
    #elif (tp=='sparse matrix'):
    #    m=re.match(r"# nnz: ([0-9]+)", fileobj.readline()); assert(m); nnz=int(m.group(1))
    #    m=re.match(r"# rows: ([0-9]+)", fileobj.readline()); assert(m); rows=int(m.group(1))
    #    m=re.match(r"# columns: ([0-9]+)", fileobj.readline()); assert(m); columns=int(m.group(1))
    #    data=read_sparse(fileobj, nnz, rows, columns)
    else: raise Exception("Unknown data type")
    assert(re.match(r"\s*", fileobj.readline()))
    assert(re.match(r"\s*", fileobj.readline()))
    return (name, data)

def read_oct(fileobj):
    if isinstance(fileobj, str):
        fileobj=open(fileobj, 'r')
    m=re.match(r"# Created by (.*)", fileobj.readline())
    res={'created':m.group(1)}
    while True:
        d=read_array(fileobj)
        if not d: break
        res[d[0]]=d[1]
    fileobj.close()
    #print()
    return res

def load_results(filename):
    octdata=read_oct(filename)
    sys=LatticeFromOct(octdata)
    path=np.transpose(octdata['PATH'], axes=(0,2,3,4,1,5))
    if path.shape[3]==1: path=path.squeeze(axis=3)
    ergy=octdata['CONTRIBUTIONS'] if 'CONTRIBUTIONS' in octdata else octdata['ENERGY']
    distance=octdata['DISTANCE'].squeeze() if 'DISTANCE' in octdata else None
    if 'GRAD' in octdata:
        grad=np.transpose(octdata['GRAD'], axes=(0,2,3,4,1,5))
        if grad.shape[3]==1: grad=grad.squeeze(axis=3)
    else: grad=None
    return (sys, path, ergy, distance, grad)

def write_ndarray(fileobj, array, name):
    if array.ndim==1: array=array[:,None]
    fileobj.write("# name: {}\n".format(name))
    fileobj.write("# type: matrix\n")
    fileobj.write("# ndims: {}\n".format(array.ndim))
    for n in array.shape: fileobj.write("{} ".format(n))
    fileobj.write("\n")
    for n in array.flatten(order='F'): fileobj.write("{}\n".format(n))
    fileobj.write("\n\n")

def write_scalar(fileobj, scalar, name):
    fileobj.write("# name: {}\n".format(name))
    fileobj.write("# type: scalar\n")
    fileobj.write("{}\n\n\n".format(scalar))

def write_lattice(fileobj, sys):
    if sys.dim==2:
        size=np.array(sys.size+(1,))
        cell=np.zeros((2,3),dtype=nptype); cell[:,:2]=np.asarray(sys.cell)
        translations=np.eye(3,dtype=nptype); translations[:2,:2]=sys.translations
        bonds=np.zeros((len(sys.bonds),5))
        bonds[:,:2]=np.array([n[0] for n in sys.bonds])
        bonds[:,3]=np.array([n[1] for n in sys.bonds])
        bonds[:,4]=np.array([n[2] for n in sys.bonds])
    else:
        size=np.array(sys.size)
        cell=np.asarray(sys.cell)
        translations=sys.translations
        bonds=np.zeros((len(sys.bonds),5))
        bonds[:,:3]=np.array([n[0] for n in sys.bonds])
        bonds[:,3]=np.array([n[1] for n in sys.bonds])
        bonds[:,4]=np.array([n[2] for n in sys.bonds])
    write_ndarray(fileobj, np.asarray(size), 'SZ')
    write_ndarray(fileobj, cell, 'CELL')
    write_ndarray(fileobj, translations, 'TRANSLATIONS')
    write_ndarray(fileobj, bonds, 'BONDS')
    write_ndarray(fileobj, sys.H, 'H')
    write_ndarray(fileobj, sys.K0, 'K0')
    write_ndarray(fileobj, sys.K, 'K')
    write_ndarray(fileobj, sys.J, 'J')
    write_ndarray(fileobj, sys.D, 'D')
    if sys.mu!=0: write_scalar(fileobj, sys.mu, 'mu')
    write_scalar(fileobj, sys.gamma, 'gamma')
    write_ndarray(fileobj, sys.bc, 'BC')

def save_results(fileobj, sys, path, energy=None, distance=None, hessian=None, grad=None):
    if isinstance(fileobj, str): fileobj=open(fileobj, 'w')
    fileobj.write("# Created by Miser library\n")
    write_lattice(fileobj, sys)
    if energy is None:
        energy_contributions=sys.lambda_energy_contributions()
        energy=np.array([energy_contributions(path[n]) for n in range(path.shape[0])])
    if distance is None:
        distance=np.cumsum([0]+[spher_distance(a,b,flat=False) for a,b in zip(path[1:],path[:-1])])
    if path.ndim==5: path=path[...,None,:,:]
    path=np.transpose(path, axes=(0,4,1,2,3,5))
    write_ndarray(fileobj, path, "PATH")
    if energy.ndim==1: write_ndarray(fileobj, energy, "ENERGY")
    else: write_ndarray(fileobj, energy, "CONTRIBUTIONS") 
    if not distance is None:
        write_ndarray(fileobj, distance, 'DISTANCE') 
    if not grad is None:
        if grad.ndim==5: grad=grad[...,None,:,:]
        grad=np.transpose(grad, axes=(0,4,1,2,3,5))
        write_ndarray(fileobj, grad, "GRAD")    
    fileobj.close()

# Graphics

def plot_against(data):
    N=data.shape[1]
    plt.subplots_adjust(hspace=0., wspace=0.)    
    ax=[]; xticklabels=[];  yticklabels=[];
    for n in range(N): 
        a=plt.subplot(N,N,1+n+n*N)
        ax.append(a)
        plt.plot(data[:,n], data[:,n])
        if n<=N-1: xticklabels.append(a.get_xticklabels())
        if n>n: yticklabels.append(a.get_yticklabels())
    for n in range(N): 
        for m in range(N): 
            if n>m:
                a=plt.subplot(N,N,1+n+m*N,sharex=ax[n],sharey=ax[m])
                a.plot(data[:,n],data[:,m])
                if m<=N-1: xticklabels.append(a.get_xticklabels())
                if n>m: yticklabels.append(a.get_yticklabels())
    plt.setp(xticklabels, visible=False)
    plt.setp(yticklabels, visible=False)

def load_dataset(folder, suffix='/mep.oct'):
    systems=[]; paths=[]; energies=[]; distances=[]; grads=[]
    folders=next(os.walk(folder))[1]
    folders.sort()
    print(folders)
    for dirname in folders:
        sys, path, ergy, distance, grad=load_results(folder+'/'+dirname+suffix)
        systems.append(sys) 
        paths.append(path) 
        energies.append(ergy)
        distances.append(distance)
        grads.append(grad)
    return systems, paths, energies, distances, grads

def save_dataset(directory, systems, paths):
    if not os.path.exists(directory): os.makedirs(directory)
    for n, (sys, path) in enumerate(zip(systems, paths)):
        dirname="{}/{}".format(directory, n)
        if not os.path.exists(dirname): os.makedirs(dirname)
        save_results("{}/mep.oct".format(dirname), sys, path)