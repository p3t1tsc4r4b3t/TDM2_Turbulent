import matplotlib.pyplot as plt
import math as m
import numpy as np
###########################################
###########################################lecture fichier
###########################################
f = open("data.txt","r")
depickler=[]
for i in f:
    depickler.append(i)
f.close()
t=[]
u=[]
v=[]
w=[]
temp=[]
for i in depickler:
    temp=i.split()
    t.append(float(temp[1]))
    u.append(float(temp[2]))
    v.append(float(temp[3]))
    w.append(float(temp[4]))

plt.figure()
plt.plot(t,u,label='u')
plt.plot(t,v,label='v')
plt.plot(t,w,label='w')
plt.xlabel("t")
plt.ylabel("vitesses")
plt.grid()
plt.legend()
#plt.savefig("uvx_t")
plt.show()
###########################################
###########################################umean
###########################################
def mean(L):
    res = 0
    n=len(L)
    for i in L:
        res+=i
    res=res/n
    return res
umean=mean(u)
vmean=mean(v)
wmean=mean(w)
###########################################
###########################################ecart type
###########################################
uprim=[]
vprim=[]
wprim=[]
n=len(u)
for i in range(n):
    uprim.append(u[i]-umean)
    vprim.append(v[i]-vmean)
    wprim.append(w[i]-wmean)
print("moyenne mesurée de u,v,w = ",umean,vmean,wmean)

uprimsquare=[]
vprimsquare=[]
wprimsquare=[]

for i in range(n):
    uprimsquare.append(uprim[i]*uprim[i])
    vprimsquare.append(vprim[i]*vprim[i])
    wprimsquare.append(wprim[i]*wprim[i])

Muprimsquare = mean(uprimsquare)
Mvprimsquare = mean(vprimsquare)
Mwprimsquare = mean(wprimsquare)
print("ecart type merusé de u,v,w : ",m.sqrt(Muprimsquare),m.sqrt(Mvprimsquare),m.sqrt(Mwprimsquare))

###########################################
###########################################tracés de densités
###########################################


###########################################
###########################################u
###########################################
xu=[0.5+0.7*i/100 for i in range(100)]
yu=[1/m.sqrt(Muprimsquare)/m.sqrt(2*m.pi)*m.exp(-0.5*((i-umean)/m.sqrt(Muprimsquare))**2) for i in xu]

plt.figure()
au,bu,su = plt.hist(u,bins=1000,density=1,label="mesures")
plt.plot(xu,yu,"r",label="courbe modèle")
plt.xlabel("u")
plt.ylabel("densité de probabilité P")
plt.grid()
plt.legend()
plt.savefig("udensity")
###########################################
###########################################v
###########################################
xv=[-0.2+0.4*i/100 for i in range(100)]
yv=[1/m.sqrt(Mvprimsquare)/m.sqrt(2*m.pi)*m.exp(-0.5*((i-vmean)/m.sqrt(Mvprimsquare))**2) for i in xv]

plt.figure()
av,bv,sv =plt.hist(v,bins=1000,density=1,stacked=True,label="mesures")
plt.plot(xv,yv,"r",label="courbe modèle")
plt.xlabel("v")
plt.ylabel("densité de probabilité P")
plt.grid()
plt.legend()
plt.savefig("vdensity")
###########################################
###########################################w
###########################################
xw=[-0.2+0.4*i/100 for i in range(100)]
yw=[1/m.sqrt(Mwprimsquare)/m.sqrt(2*m.pi)*m.exp(-0.5*((i-wmean)/m.sqrt(Mwprimsquare))**2) for i in xw]

plt.figure()
aw,bw,sw = plt.hist(w,bins=1000,density=1,stacked=True,label="mesures")
plt.plot(xw,yw,"r",label="courbe modèle")
plt.xlabel("w")
plt.ylabel("densité de probabilité P")
plt.grid()
plt.legend()
plt.savefig("wdensity")

def moycalc(a,b):
    res=0
    l=len(a)
    for i in range(l):
        res = res + a[i]*(b[i+1]-b[i])*(b[i+1]+b[i])/2
    return res

def varcalc(a,b):
    res=0
    l=len(a)
    for i in range(l):
        res = res + a[i]*(b[i+1]-b[i])*(b[i+1]+b[i])**2/4
    return res - moycalc(a,b)**2

print("moyenne recalculée de u, v, w : ", moycalc(au,bu),moycalc(av,bv),moycalc(aw,bw))
print("ecart type recalculé de u, v, w : ", m.sqrt(varcalc(au,bu)),m.sqrt(varcalc(av,bv)),m.sqrt(varcalc(aw,bw)))


###########################################
###########################################dérivée uvw
###########################################
dt=t[2]-t[1]
tprim=[dt/2+t[i] for i in range(n-1)]
uder=[(u[i+1]-u[i])/dt for i in range(n-1)]
vder=[(v[i+1]-v[i])/dt for i in range(n-1)]
wder=[(w[i+1]-w[i])/dt for i in range(n-1)]
###########################################
plt.figure()
plt.plot(tprim,uder,label='du/dt')

plt.xlabel("t")
plt.ylabel("dérivée temporelle")
plt.grid()
plt.legend()
plt.savefig("du_t")

###########################################
plt.figure()

plt.plot(tprim,vder,label='dv/dt')

plt.xlabel("t")
plt.ylabel("dérivée temporelle")
plt.grid()
plt.legend()
plt.savefig("dvdt_t")

###########################################
plt.figure()

plt.plot(tprim,wder,label='dw/dt')
plt.xlabel("t")
plt.ylabel("dérivée temporelle")
plt.grid()
plt.legend()
plt.savefig("dwdt_t")

###########################################
###########################################du mean
###########################################
print("moyenne mesurée de dudt, dvdt, dwdt : ", mean(uder),mean(vder),mean(wder))
###########################################
########################################### densité
###########################################


###########################################
###########################################du densité
###########################################

plt.figure()
plt.hist(uder,bins=1000,density=1,stacked=True,label="mesures")
plt.xlabel("du/dt")
plt.ylabel("densité de probabilité P")
plt.grid()
plt.legend()
plt.savefig("dudtdensity")
###########################################
###########################################dv densité
###########################################
plt.figure()
plt.hist(vder,bins=1000,density=1,stacked=True,label="mesures")
plt.xlabel("dv/dt")
plt.ylabel("densité de probabilité P")
plt.grid()
plt.legend()
plt.savefig("dvdtdensity")
###########################################
###########################################dw densité
###########################################
plt.figure()
plt.hist(wder,bins=1000,density=1,stacked=True,label="mesures")
plt.xlabel("dw/dt")
plt.ylabel("densité de probabilité P")
plt.grid()
plt.legend()
plt.savefig("dwdtdensity")


