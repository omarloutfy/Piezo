import sympy
from sympy import Array 
from sympy import *
import numpy as np
import sys

    
    
R=3 
    
class Tensor :
    
    
    
    def __init__(self, n, d):
            
        self.n = n
        self.d = d
        # print('HI')
    
    def Tens(self,n,d):     #Methode qui prend une liste des composante et qui retourne un tenseur d'ordre n et de dimension d
    
        def rec1(n,d):      #Methode qui prend n et d et retourne une liste des composantes d'un tenseur d'ordre n et dimension d
            
            q='P'
            u=''
            z=[]
            c=[]
            f=[] 
            k=0
            s=0
            
            if n==1:
                
                a=[]         
                for i in range(n):
                    a.append(d)
                    
                for i in range(d):
                    z=[i+1]
                    a='P{}'.format(*z)
                    c.append(sympy.core.symbol.Symbol(a))
                P=MutableSparseNDimArray(c,(d,1))
                return c
                          
                
            else:
                

                for i in range(d**(n-1)):
                    for j in range(d):
                        
                        x=str(rec1(n-1,d)[i])+'{}'.format(j+1)   
                        f.append(sympy.core.symbol.Symbol(x))           
                return f
            
            
        a=[]
        for i in range(n):
            a.append(d)   
        if n==1:
            P=MutableSparseNDimArray(rec1(1,d),[d,1])
        else :
            P=MutableSparseNDimArray(rec1(n,d),a)
        
        return P
    
    # def Input_Tensor(self):    #Tenseur ordre n et dim d
        
    #     r=self.rec1(self.n,self.d)
    #     a=[]
        
    #     for i in range(self.n):
    #         a.append(self.d)
        
    #     P=MutableSparseNDimArray(r,a)
    #     # print(P)
    #     return P
    
    def Tensor3_Sym(self):   #Tenseur ordre 3 totalement symetrique
        
        R=3
        
        P=self.Tens(3,3)
        for i in range (R):
            for j in range (R):
                for k in range (R):
        
                    P[i,j,k]=P[i,k,j]=P[k,j,i]=P[j,i,k]
        return P
    
    def Piezo_Tensor3(self):   # Tenseur piezo symetrique % 2 derniers indices
        
        P=self.Tens(3,3)
        
        for i in range (R):
            for j in range (R):
                for k in range (R):
        
                    P[i,j,k]=P[i,k,j]
        
        # print('\n[P] =', P)
        return P
    
    def Tens3_To_MatKelv_P(self):     #Tenseur piezo to matrice de Kelvin
        
        P = self.Piezo_Tensor3()
        c=([P[0,0,0],P[0,0,1],P[0,2,2],sqrt(2)*P[0,1,2],sqrt(2)*P[0,0,2],sqrt(2)*P[0,0,1],[P[1,0,0],P[1,1,1],P[1,2,2],sqrt(2)*P[1,1,2],sqrt(2)*P[1,0,2],sqrt(2)*P[1,0,1]],[P[2,0,0],P[2,1,1],P[2,2,2],sqrt(2)*P[2,1,2],sqrt(2)*P[2,0,2],sqrt(2)*P[2,0,1]]])
        p = MutableSparseNDimArray(c,(3,6)) 
        
        # print('\n[p] =', p)
        
        return p 
    
    def Tens3_To_MatKelv(self, P):    #Tenseur ordre 3 quelconques to matrice de Kelvin
            
            
        c=([P[0,0,0],P[0,0,1],P[0,2,2],sqrt(2)*P[0,1,2],sqrt(2)*P[0,0,2],sqrt(2)*P[0,0,1],[P[1,0,0],P[1,1,1],P[1,2,2],sqrt(2)*P[1,1,2],sqrt(2)*P[1,0,2],sqrt(2)*P[1,0,1]],[P[2,0,0],P[2,1,1],P[2,2,2],sqrt(2)*P[2,1,2],sqrt(2)*P[2,0,2],sqrt(2)*P[2,0,1]]])
        p = MutableSparseNDimArray(c,(3,6)) 
            
        # print('\n[p] =', p)
            
        return p 
        
        
    def Tens3_To_Voigt(self, P):    #Tenseur piezo to matrice de Voigt
        
       
        
        d = ([P[0,0,0],P[0,0,1],P[0,2,2],P[0,1,2],P[0,0,2],P[0,0,1],[P[1,0,0],P[1,1,1],P[1,2,2],P[1,1,2],P[1,0,2],P[1,0,1]],[P[2,0,0],P[2,1,1],P[2,2,2],P[2,1,2],P[2,0,2],P[2,0,1]]])
        pv = MutableSparseNDimArray(d,(3,6)) 
        # print('\n[pv] =', pv)
        return pv
    
    
    def Tens3_To_Voigt_P(self):    #Tenseur ordre 3 quelconques to matrice de Voigt
        
        P = self.Piezo_Tensor3()
        d = ([P[0,0,0],P[0,0,1],P[0,2,2],P[0,1,2],P[0,0,2],P[0,0,1],[P[1,0,0],P[1,1,1],P[1,2,2],P[1,1,2],P[1,0,2],P[1,0,1]],[P[2,0,0],P[2,1,1],P[2,2,2],P[2,1,2],P[2,0,2],P[2,0,1]]])
        pv = MutableSparseNDimArray(d,(3,6)) 
        # print('\n[pv] =', pv)
        return pv
        
class Transform_CG(Tensor) :   #Class pour la decomposition de Clebsh-Gordan
            
    def __init__(self):
        super().__init__(3,3)
        
####################################################################################       
############ Methodes qui prennent P et retournent les harmoniques ################

    def Harm_a (self) :      
        
        
        P = self.Piezo_Tensor3()
        l2=[]
        for i in range(R):  
            
            for j in range(R) :
                z1=0
                for p in range(R):
                    
                    for q in range(R):
                        
                        z1=factor(z1+(Rational(1,2)*(LeviCivita(i,p,q)*P[p,q,j]+(LeviCivita(j,p,q)*P[p,q,i]))))                
                l2.append(z1)
        a = MutableSparseNDimArray(l2, (R,R))  
        
        # print('\n[a] =', a)
        return a      
    
    
    def Harm_u (self) :
            
            P = self.Piezo_Tensor3()
            l3=[]
            for i in range(R): 
                z2=0
                for p in range(R):
                   
                    z2=factor(z2+P[p,p,i]-(Rational(1,3))*P[i,p,p])                
                l3.append(z2)                
            u = Array(l3, (R,1))  
            
            # print('\n[u] =', u)
            return u
        
        
    def Harm_v (self) :
        
        P = self.Piezo_Tensor3()
        l4=[]
        for i in range(R): 
            z3=0
            for p in range(R):
               
                z3=z3+P[i,p,p]
            
            l4.append(z3)
            
        v = Array(l4, (R,1))
        # print('\n[v] =', v)
        return v
        
        
        
    def Harm_H (self) :
        
        P = self.Piezo_Tensor3()
        a = self.Harm_a()
        u = self.Harm_u()
        v = self.Harm_v()
        l5=[]
        l6=[]
        
        for i in range(R):  
             
            for j in range(R) :
                
                for k in range(R):
                    z5=0
                    z4=0
                    z6=0
                    z7=0
                    for p in range(R):        
                        
                        z4=(z4-Rational(1,3)*(LeviCivita(i,j,p)*a[p,k]+LeviCivita(i,k,p)*a[p,j]))
                        z6=(-Rational(3,10)*(u[j,0]*KroneckerDelta(i,k)+u[k,0]*KroneckerDelta(i,j)-Rational(2,3)*u[i,0]*KroneckerDelta(j,k)))
                        z7=(-Rational(1,3)*v[i,0]*KroneckerDelta(j,k)   )
                        
                        z5=factor(P[i,j,k]+z4+z6+z7)
                          
                    l5.append(z5)
                    l6.append(z4)
                                                   
        H = Array(l5, (R,R,R))   #affichage matrice a
        # print('\n[H] =', H)
        return H

    

    def Harm_P12(self) :
        
        H = self.Harm_H()
        a = self.Harm_a()
        u = self.Harm_u()
        l7=[]
        
        for i in range(R):  
             
            for j in range(R) :
                
                for k in range(R):
                    z5=0
                    z4=0
                    z6=0
                    z7=0
                    for p in range(R):        
                        
                        # r=a.append(((0.5)*(LC[i,p,q]*P[p,q,j]+(LC[j,p,q]*P[p,q,i])))) 
                        z4=z4+Rational(1,3)*(LeviCivita(i,j,p)*a[p,k]+LeviCivita(i,k,p)*a[p,j])
                        z6=Rational(3,10)*(u[j,0]*KroneckerDelta(i,k)+u[k,0]*KroneckerDelta(i,j)-Rational(2,3)*u[i,0]*KroneckerDelta(j,k))
                         
                        
                        z8=factor(H[i,j,k]+z4+z6)
            
                    l7.append(z8)
        
        P12 = Array(l7, (R,R,R)) 
        # print('\n[P12] =', P12)
        return P12
    
    
    def Harm_P10(self) :
         
        v = self.Harm_v()
        l8=[]
        
        for i in range(R):  
             
            for j in range(R) :
                
                for k in range(R):
                    z5=0
                    z4=0
                    z6=0
                    z7=0
                    for p in range(R):        
                        
                        # r=a.append(((0.5)*(LC[i,p,q]*P[p,q,j]+(LC[j,p,q]*P[p,q,i])))) 
                        z7=Rational(1,3)*v[i,0]*KroneckerDelta(j,k)                                         
                        z9=factor(z7)  
                   
                    l8.append(z9)
        
        P10 = Array(l8, (R,R,R))   
        # print('\n[P10] =', P10)
        return P10

#####################################################################################
### Methodes qui prend les harmonique et normallement retourne le tenseur Piezo P ###


    def P_CG(self) :
        
        H = self.Harm_H()
        a = self.Harm_a()
        u = self.Harm_u()
        v = self.Harm_v()
        
        n = MutableSparseNDimArray(range(50),(R,R,R))  #initialisations
        r = MutableSparseNDimArray(range(50),(R,R,R))
        w = MutableSparseNDimArray(range(50),(R,R,R))
        
        n=simplify(P12+P10)
        
        
        
        l8=[]
        for i in range (R):    
            for j in range (R):
                for k in range (R):
                    z1=0
                    for q in range(R):
                        z1 = z1+LeviCivita(i,j,q)*a[q,k]+LeviCivita(i,k,q)*a[q,j]
                    r[i,j,k] = H[i,j,k]+Rational(1,3)*z1+Rational(3,10)*((u[j,0])*KroneckerDelta(i,k)+u[k,0]*KroneckerDelta(i,j)-Rational(2,3)*u[i,0]*KroneckerDelta(j,k))+Rational(1,3)*v[i,0]*KroneckerDelta(j,k)
                    s=simplify(r)     
                    
                    
                    # w[i,j,k]=s[i,j,k]-P[i,j,k]   
                    
        return s  

    def Tr3(self,H): #Trace d'un tenseur ordre 3 (H[k,k,i]) placé dans un vecteur
        
        a=[]
        for i in range (R):
            ss=0
            for j in range(R):
                ss=simplify(H[j,j,i])+ss
            a.append(ss)
            v=MutableSparseNDimArray(a,(3,1))
        return v
            
        



class Transform_SW(Tensor) : #Class pour la decomposition de Schur-Weyl
            
    def __init__(self):
        super().__init__(3,3)
        

####################################################################################       
############ Methodes qui prennent P et retournent les harmoniques ################

    def Harm_a (self) :
        
        
        P = self.Piezo_Tensor3()
        l2=[]
        for i in range(R):  
            
            for j in range(R) :
                z1=0
                for p in range(R):
                    
                    for q in range(R):
                        
                        z1=factor(z1+(Rational(1,2)*(LeviCivita(i,p,q)*P[p,q,j]+(LeviCivita(j,p,q)*P[p,q,i]))))                
                l2.append(z1)
        a = MutableSparseNDimArray(l2, (R,R))  
        
        # print('\n[a] =', a)
        return a      
    
    
    def Harm_u (self) :
            
            P = self.Piezo_Tensor3()
            l1=[]
            for i in range(R): 
                z1=0
                for p in range(R):
                   
                    z1=z1+Rational(1,3)*(P[i,p,p]+2*P[p,p,i])
                
                l1.append(z1)
                
            u = Array(l1, (R,1))  
            
            # print('\n[u] =', u)
            return u
        
        
    def Harm_v (self) :
        
        P = self.Piezo_Tensor3()
        l2=[]
        for i in range(R): 
            z2=0
            for p in range(R):
               
                z2=z2+Rational(1,3)*(P[p,p,i]-P[i,p,p])
            
            l2.append(z2)
            
        v = Array(l2, (R,1))  
        # print('\n[v] =', v)
        return v
        
        
        
    def Harm_H (self) :
        
        P = self.Piezo_Tensor3()
        a = self.Harm_a()
        u = self.Harm_u()
        v = self.Harm_v()
        l4=[]
        z4=0
        for i in range(R):  
             
            for j in range(R) :
                
                for k in range(R):
                    
                    z4=Rational(1,3)*(P[i,j,k]+P[k,i,j]+P[j,i,k])-Rational(1,5)*(KroneckerDelta(i,j)*u[k,0]+KroneckerDelta(i,k)*u[j,0]+KroneckerDelta(j,k)*u[i,0])
                    l4.append(z4)
        
        
        H = Array(l4,(R,R,R))
        # print('\n[H] =', H)
        return H




    def Harm_Ps(self) :
        
        H = self.Harm_H()
        u = self.Harm_u()
        
        l5=[]
        z5=0
        for i in range(R):  
             
            for j in range(R) :
                
                for k in range(R):
                    
                    z5=H[i,j,k]+Rational(1,5)*(KroneckerDelta(i,j)*u[k,0]+KroneckerDelta(i,k)*u[j,0]+KroneckerDelta(j,k)*u[i,0])
                    l5.append(z5)
        
        Ps = Array(l5,(R,R,R))
        # print('\n[Ps] =', Ps)
        return Ps
    
    
    def Harm_Pr(self) :
         
        a = self.Harm_a()
        v = self.Harm_v()
        l8=[]
        for i in range(R):  
             
            for j in range(R) :
                
                for k in range(R):
        
                    z6=0
                    for p in range(R):        
                        
                        # r=a.append(((0.5)*(LC[i,p,q]*P[p,q,j]+(LC[j,p,q]*P[p,q,i])))) 
                        z6=z6+Rational(1,3)*(LeviCivita(i,j,p)*a[p,k]+LeviCivita(i,k,p)*a[p,j])
                        z7=-Rational(1,2)*(2*v[i,0]*KroneckerDelta(j,k)-v[k,0]*KroneckerDelta(i,j)-v[j,0]*KroneckerDelta(i,k))
                        
                        z8=z6+z7
                        
                    #print(z4)   
                    #print(z5)   
                    
                    l8.append(z8)
        

        Pr = Array(l8,(R,R,R))
    
        # print('\n[Pr] =', Pr)
        return Pr
    
#####################################################################################
### Methodes qui prend les harmonique et normallement retourne le tenseur Piezo P ###
    
                
    def P_SW(self) : #construction de P a partir des harmoniques H a 
        
        H = self.Harm_H()
        a = self.Harm_a()
        u = self.Harm_u()
        v = self.Harm_v()
        
        o = MutableSparseNDimArray(range(50),(R,R,R))  #initialisations
        l = MutableSparseNDimArray(range(50),(R,R,R))
        t = MutableSparseNDimArray(range(50),(R,R,R))
        
        o=simplify(Ps+Pr)
        
        
        
        l8=[]
        for i in range (R):    
            for j in range (R):
                for k in range (R):
                    z1=0
                    for q in range(R):
                        z1 = z1+LeviCivita(i,j,q)*a[q,k]+LeviCivita(i,k,q)*a[q,j]
                    l[i,j,k] = H[i,j,k]+Rational(1,3)*z1+Rational(1,5)*((u[j,0])*KroneckerDelta(i,k)+u[k,0]*KroneckerDelta(i,j)+u[i,0]*KroneckerDelta(j,k))-Rational(1,2)*(2*v[i,0]*KroneckerDelta(j,k)-v[k,0]*KroneckerDelta(i,j)-v[j,0]*KroneckerDelta(i,k))
                    c=simplify(l)     
                    
                    
                    t[i,j,k]=c[i,j,k]-P[i,j,k]   

        return c
    
    
    def Tr3(self,H): #Trace d'un tenseur ordre 3 (H[k,k,i]) placé dans un vecteur
        
        a=[]
        for i in range (R):
            ss=0
            for j in range(R):
                ss=simplify(H[j,j,i])+ss
            a.append(ss)
            v=MutableSparseNDimArray(a,(3,1))
        return v
    
###################################################################################


r = Tensor(10,10)

m = Tensor(3,3)

P1 = m.Tens(1,4)    #Tenseur ordre 1 dim 4
P2 = m.Tens(2,8)    #Tenseur ordre 2 dim 8
P3 = m.Tens(3,2)    #Tenseur ordre 3 dim 2
P4 = m.Tens(4,3)    #Tenseur ordre 4 dim 3
P5 = m.Tens(5,2)    #Tenseur ordre 5 dim 2

T3 = m.Tensor3_Sym()      #Tenseur ordre 3 Totalement symetrique
Tp = m.Piezo_Tensor3()   #Tenseur piezo





P = r.Piezo_Tensor3()            #Tenseur Piezo

p = r.Tens3_To_MatKelv_P()       #Mat kelvin p tiré du tenseu P piezo
pv = r.Tens3_To_Voigt_P()        #Mat voigt p tiré du tenseu P piezo

p1 = r.Tens3_To_MatKelv(P)       #Mat kelvin tiré de tenseur ordre 3 quelconques P
p2 = r.Tens3_To_Voigt(P)         #Mat voigt tiré de tenseur ordre 3 quelconques P




####################################################################################
    

################ Harmoniques avec la décomposition de Clebsh-Gordan ################

d1 = Transform_CG()  
     
a = d1.Harm_a()
# u = d1.Harm_u()
# v = d1.Harm_v()
H = d1.Harm_H()  

TrCG=d1.Tr3(H)      #Trace de l'harmonique H sous Clebsh-Gordan placé dans un vecteur


# KelvH = d1.Tens3_To_MatKelv(H)   #test Tens3_To_MatKelv d'un tenseur 3 qqoncque 
VoiH = d1.Tens3_To_Voigt(H)

P10 = d1.Harm_P10()
P12 = d1.Harm_P12()

PCG = d1.P_CG()





######### consrtuction de P a partir des harmoniques de Clebsh-Gordan #########


# n = MutableSparseNDimArray(range(50),(R,R,R))  #initialisations
# r = MutableSparseNDimArray(range(50),(R,R,R))
# w = MutableSparseNDimArray(range(50),(R,R,R))

# n=simplify(P12+P10)



# l8=[]
# for i in range (R):    
#     for j in range (R):
#         for k in range (R):
#             z1=0
#             for q in range(R):
#                 z1 = z1+LeviCivita(i,j,q)*a[q,k]+LeviCivita(i,k,q)*a[q,j]
#             r[i,j,k] = H[i,j,k]+Rational(1,3)*z1+Rational(3,10)*((u[j,0])*KroneckerDelta(i,k)+u[k,0]*KroneckerDelta(i,j)-Rational(2,3)*u[i,0]*KroneckerDelta(j,k))+Rational(1,3)*v[i,0]*KroneckerDelta(j,k)
#             s=simplify(r)     #s = tenseur 3 P normallement
            
            
#             w[i,j,k]=s[i,j,k]-P[i,j,k]   #tenseur nulle normallement
            

# print('\n[n] =', n)     # P12+P10=P                               vérifié
# print('\n[s] =', s)     # Decomp harmonique forme le tenseur P    vérifié


#ectriture avec Rational(a,b) permet d'eviter le 0 numérique

#Resultat : on retrouve bien le tenseur P donc la decomposition a été construite correctement
# dans le cas de Clebsch-Gordan

####################################################################################


############## Harmoniques avec la décomposition de Schur-Weyl ##############

d2 = Transform_SW()    

# aa = d2.Harm_a()
# uu = d2.Harm_u()
# vv = d2.Harm_v()
HH = d2.Harm_H()  

TrSW=d2.Tr3(HH)      #Trace de l'harmonique H sous Clebsh-Gordan placé dans un vecteur


Ps = d2.Harm_Ps()
Pr = d2.Harm_Pr()
PSW = d2.P_SW()

########## consrtuction de P a partir des harmoniques de Schur-Weyl  ##########


# o = MutableSparseNDimArray(range(50),(R,R,R))  #initialisations
# l = MutableSparseNDimArray(range(50),(R,R,R))
# t = MutableSparseNDimArray(range(50),(R,R,R))

# o=simplify(Ps+Pr)



# l8=[]
# for i in range (R):    
#     for j in range (R):
#         for k in range (R):
#             z1=0
#             for q in range(R):
#                 z1 = z1+LeviCivita(i,j,q)*a[q,k]+LeviCivita(i,k,q)*a[q,j]
#             l[i,j,k] = HH[i,j,k]+Rational(1,3)*z1+Rational(1,5)*((uu[j,0])*KroneckerDelta(i,k)+uu[k,0]*KroneckerDelta(i,j)+uu[i,0]*KroneckerDelta(j,k))-Rational(1,2)*(2*vv[i,0]*KroneckerDelta(j,k)-vv[k,0]*KroneckerDelta(i,j)-vv[j,0]*KroneckerDelta(i,k))
#             c=simplify(l)     #c = tenseur 3 P normallement
            
            
#             t[i,j,k]=c[i,j,k]-P[i,j,k]   #tenseur nulle normallement si decomp correcte
            
            
# print('\n[o] =', np.array(o))     # Ps+Pr=P                                   vérifié
# print('\n[c] =', np.array(c))     # Decomp harmonique forme le tenseur P      vérifié



#Resultat : on retrouve bien le tenseur P donc la decomposition a été construite correctement
# dans le cas de Schur-Weyl



# Sur console, taper  simplify(P10+P12)   ou    s pour decomp CG
# Sur console, taper  simplify(Ps+Pr)     ou    c pour decomp SW

s=0
for i in range(R):
        for j in range(R):
            for k in range(R):
                s=s+a[k,k]
            if (a[i,j]==a[j,i]) and s==0:
                a[1,2]=a[2,1]=sympy.core.symbol.Symbol("a{}".format(3))
                a[0,2]=a[2,0]=sympy.core.symbol.Symbol("a{}".format(4))
                a[0,1]=a[1,0]=sympy.core.symbol.Symbol("a{}".format(5))
                for l in range(R):
                    a[l,l]=sympy.core.symbol.Symbol("a{}".format(l+1))
                    a[2,2]=-a[0,0]-a[1,1]
                
    
       
        
            

