
  
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 15:59:36 2019
@author: tomsa
"""
#%%
import numpy as np 
import copy
import numpy.random as rnd
import copy
import numpy.random as rnd
import time
import matplotlib.pyplot as plt
import pandas as pd


#%%
def cout_soudure(x):
    return 1.7781*0.625*x[1]*x[2]**2 + 0.6224*0.625*x[0]*x[2]*x[3] + 3.1661*(0.625*x[0])**2*x[3] + 19.84*(0.625*x[0])**2*x[2]

#%%
def g(i,x): #contraintes
    return { 1 : 0.00954*x[2]-0.625*x[1], \
             2 :  0.0193*x[2]- 0.625*x[0], \
             3 : -np.pi * x[2]**2 * x[3] - 4/3 * np.pi * x[2]**3 + 1296000 }[i]

#%%
def penalite(x):
    p = [max(g(i,x),0) for i in range(1,4)]

    return np.sum(p)


#%%
def initiation(f,bounds,p):
   

    d=len(bounds) #finding number of dimensions
    pop=np.zeros(p) #creating empty position array
    pop=pop.tolist() #converting array to list
    pop_velocity=pop[:] #empty velocity array
    pop_val=pop[:] #empty value array
    for j in range(p): #iterating ovre the number of particles
        pop[j]=[rnd.uniform(bounds[i][0],bounds[i][1])\
                    for i in range(d)] #random coordinate within bounds
        pop_val[j]=f(pop[j]) #calculating function value
                                            #at each particle
        pop_velocity[j]=[rnd.uniform(-abs(bounds[i][1]-bounds[i][0])\
                    ,abs(bounds[i][1]-bounds[i][0])) for i in range(d)]
                    #creating random velocity values for each dimension

    global_best=pop[np.argmin(pop_val)]#getting the lowest particle value
    pop_best=copy.deepcopy(pop)#setting all particles current positions to best
    return d,np.array(pop), np.array(pop_best), \
                 np.array(global_best), np.array(pop_velocity), \
                     np.array(pop_val)


#%%       
def withinbounds(bounds,pop):
   
    for i in range(len(bounds)):
        if pop[i]<bounds[i][0]: #if particle is less than lower bound
            pop[i]=bounds[i][0]
        elif pop[i]>bounds[i][1]: #if particle is more than higher bound
            pop[i]=bounds[i][1]
    return

#%%
def print_results(hist):
    
    scores = hist[:,1]
    best_score = min(scores)
    index = np.where(scores == best_score)
    #print(f"index : {index}")
    print(f"meilleurs dimensions : \n z1 = {0.625*hist[index[0][0],3]} \n z2 = {0.625*hist[index[0][0],4]} \n x3 = {hist[index[0][0],5]} \n x4 = {hist[index[0][0],6]}\n")
    print(f"meilleur score : {best_score}\n")

   
# %%
def preparation_table(hist):
    hist[:,3] = [0.625 * a for a in hist[:,3]]
    hist[:,4] = [0.625 * a for a in hist[:,4]]

    table = pd.DataFrame(data=hist[:,1:], columns=["Score","Nbr itrs", "z1","z2","x3","x4"], index=hist[:,0]) 
    
    return table




# %%
def print_table(hist):
      print("historique detaille :\n")
      return preparation_table(hist)


# %%
def print_analyse(hist):
    print("Analyse des donnees generees :\n")
    return preparation_table(hist).describe()

#%%
def particleswarm(f,exec, iter, bounds,p,c1,c2,vmax,tol):

   
    c3=c1+c2
    K=2/(abs(2-c3-np.sqrt((c3**2)-(4*c3)))) #creating velocity weighting factor
    
    historique = np.empty([0,7])
    
    count = 0
    count1=0
    for k in range (exec):

        d,pop, pop_best, global_best, pop_velocity, \
        pop_best_scores \
        = initiation(f,bounds,p) #initializing various arrays
        old_global_best=[0]*d
     

        print("global best IS :", global_best, count1, count)
        count = 0
        count1=0
        for j in range(iter): 
    
    
                if abs(f(old_global_best)-f(global_best))>tol: #exit condition
                    total_iter =j
                    break 

                               
                for i in range(p): #iterates over each particle
                    rp,rg=rnd.uniform(0,1,2) #creates two random numbers between 0-
                    pop_velocity[i,:]+=(c1*rp*(pop_best[i,:]-pop[i,:]))
                    pop_velocity[i,:]+=(c2*rg*(global_best[i,:]-pop[i,:]))
                    pop_velocity[i,:]=pop_velocity[i,:]*K
                    if pop_velocity[i].any() > vmax : #is any velocity is greater than vmax
                            pop_velocity[i,:]=vmax #set velocity to vmax
                    #all of the above is regarding updating the particle's velocity
                    #with regards to various parameters (local_best, p_best etc..)
                    pop[i,:]+=pop_velocity[i,:] #updating position
                    
                    withinbounds(bounds,pop[i]) #if particle is out of bounds

                    penal = penalite(pop[i])
                    particle_fitness= f(pop[i]) + penal # penalite = 0 si la particule respecte les contraites
                
                    
                    if particle_fitness < pop_best_scores[i]:
                        pop_best[i,:]=pop[i,:] #checking if new best
                        pop_best_scores[i]=particle_fitness
                        f_global_best=f(global_best)
                        count1+=1
                        if particle_fitness < f_global_best : 
                            count+=1
                            old_global_best=global_best[:]
                            global_best=copy.deepcopy(pop_best[i,:]) 
                            print('current function value: ',f(global_best))


    
    #print('Optimum at: ',global_best,'\n','Function at optimum: ',f(global_best)) 

        historique = np.append(historique,[np.append([k,f(global_best),total_iter if total_iter > 0 else iter],[global_best])],axis=0)
    
    return historique

#%%               
f=cout_soudure

bounds = [[1,99],[1,99],[10,200],[10,240]]    
p=60 #shouldn't really change 
vmax=(bounds[3][1]-bounds[3][0])*0.75

c1=2.8 #shouldn't really change
c2=1.3 #shouldn't really change
tol=0.0001
exec =50
iter = 100
#%%
hist = particleswarm(f,exec,iter,bounds,p,c1,c2,vmax,tol)    

# %%
print_results(hist)


# %%
print_table(hist)


# %%
print_analyse(hist)

                
                
 ###REFERENCES
 ###https://www.cs.cinvestav.mx/~constraint/papers/eisci.pdf
 ###https://github.com/TomRSavage/ParticleSwarm/blob/master
 ###https://www.sciencedirect.com/science/article/pii/S0957417419305925

# %%
