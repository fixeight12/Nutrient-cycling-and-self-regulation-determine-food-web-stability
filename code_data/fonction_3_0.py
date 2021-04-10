# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:45:01 2020

@author: Theis Kevin
"""

import numpy as np 

import csv

import os



#  0 => Non/Faux/False 1 => Oui/Vrai/True


# Change Current Working Directory

path = "C:\\Users\\Theis Kevin\\Documents\\MASTER Sciences de la Mer\\Stage M2\\projets\\david_shanafelt3.0"
os.chdir(path)


file = open("D_TL1_pertuB0.csv","r") # fichier des parametres / Parameters File

reader = csv.reader(file,delimiter = ";", skipinitialspace=True)
header = next(reader) # une fois donc premiere ligne de reader sont les headers
dictio = []
for row in reader:
    dictio.append(dict(zip(header, map(float, row))))





# dict() => creation dictionnaire / zip => permet associer un element d'une liste a une autre
# map permet d'appliquer une fonction ici float sur un element iterable ici les lignes


###############################################################################



#Création des matrices / vecteurs utilisée
    
# Différente valeurs de a (taux de predation) et b (biomass conversion efficiency)
# Value of a (attack rate) and b ()    

def vect_a_b(n,par_1): 
    a = 0 * np.ones(4)
    b = 0 * np.ones(4)
    for i in range(4):
        if i == 0 :
            a[i] = par_1["a2"]
            b[i] = par_1["b2"]
            
        elif i == 1:
            a[i] = par_1["a3"]
            b[i] = par_1["b3"]
            
        elif i == 2:
            a[i] = par_1["a4"]
            b[i] = par_1["b4"]
        
        else:
            a[i] = par_1["a5"]
            b[i] = par_1["b5"]
    return a[:n] , b[:n]
    



def vect_gi(B1,n,par_1):
    gi = 0 * np.ones(n) # creation du vecteur de croissance/mort independante de la densité
    for i in range (n):
        if i == 0:
            gi[i] =  par_1["a2"] * B1
        else:
            gi[i] = 0 
    return gi
            


def vecteur_para(n,par_1) :
    ratio_1CN = 0 * np.ones(n) # vecteur de 1/C:N
    ratio_1CN_no_convers = 0 * np.ones(n)
    for i in range (n):  
        if i == 0:
            ratio_1CN[i] = 1/par_1["alpha2"]
            ratio_1CN_no_convers[i] = 1/par_1["alpha2"]
            
        elif i == 1 :
            ratio_1CN[i] = 1 / par_1["alpha3"]
            ratio_1CN_no_convers[i] = 1/((par_1["alpha3"]*par_1["alpha2"]*(1-par_1["b3"]))/(par_1["alpha3"]-par_1["alpha2"]*par_1["b3"]))
        
        else:
            ratio_1CN[i] = 1 / par_1["alpha3"]
            ratio_1CN_no_convers[i] = 1 / par_1["alpha3"]
    return ratio_1CN , ratio_1CN_no_convers
            
        

def matrix_para(n,par_1):  
    predation = np.diag(vect_a_b(n,par_1)[1][1:]*vect_a_b(n,par_1)[0][1:]*np.ones(n-1),k = -1) # matrice de conversion de la proie // prey conversion matrix
    prey = np.diag(vect_a_b(n,par_1)[0][1:]*np.ones(n-1), k = 1) # matrice de mort // Death matrix
    no_conversion = np.diag((1-vect_a_b(n,par_1)[1][1:])*vect_a_b(n,par_1)[0][1:]*np.ones(n-1),k = -1)
    return predation , prey, no_conversion



    
    
####  Fonction de recyclage // Recycling function

def recyclage(B,n,par_1):     # un verteur de taille n doit rentrer et attention B1 = N ICI
   
    recy_mort = par_1["lam"]*vecteur_para(n,par_1)[0] * (par_1["m"] + par_1["D"] * B) * B
    recy_noconvers = par_1["lam"] * vecteur_para(n,par_1)[1] * np.dot(matrix_para(n,par_1)[2],B) * B
    
    recy_mort_tot = sum(recy_mort)
    recy_noconvers_tot = sum(recy_noconvers)
    recy_tot = recy_mort_tot + recy_noconvers_tot
    recyclage_par_sp = []
    recyclage_sp_mort =[]
    recyclage_sp_noconvers =[]
    for i in range(len(B)):
        recyclage_par_sp.append(recy_mort[i]+recy_noconvers[i])
        recyclage_sp_mort.append(recy_mort[i])
        recyclage_sp_noconvers.append(recy_noconvers[i])
        
    recy_tot = np.array([recy_tot])
    recy_mort_tot = np.array([recy_mort_tot])
    recy_noconvers_tot = np.array([recy_noconvers_tot])
    W = np.concatenate([recy_tot,recy_mort_tot,recy_noconvers_tot,recyclage_par_sp,recyclage_sp_mort,recyclage_sp_noconvers])
    return W
    

# Création du systeme // creation of the system ###
    
def syst(X,t,n,par_1):
    B1 = X[0]
    B = X[1:]
    dN = par_1["I"] - par_1["l"] * B1 - par_1["a2"] / par_1["alpha2"] * B1 * B[0] + recyclage(B,n,par_1)[0]
    Y = B * (vect_gi(B1,n,par_1) + np.dot(matrix_para(n,par_1)[0],B) - np.dot(matrix_para(n,par_1)[1],B) - par_1["m"] - par_1["D"] * B)
    dN = np.array([dN])
    W = np.concatenate([dN,Y])
    return W


def flux(X,t,n,par_1):
    B1 = X[0]
    B = X[1:]
    # Flux pour B1 // Flow for B1
    flux_I = par_1["I"]
    flux_mortB1 = par_1["l"]* B1
    flux_huntB1 = par_1["a2"]/par_1["alpha2"]*B1*B[0]
    flux_recy = recyclage(B,n,par_1)[3:(3+n)]
    
    # Flux pour autre Bi (=! B1) // Flow for Others
    flux_mortB = B * (par_1["m"] + par_1["D"]*B)
    
    flux_prelB2 = par_1["a2"]*B1*B[0]
    flux_prelB = B * np.dot(matrix_para(n,par_1)[0],B)
    flux_prelB = flux_prelB[1:]
    flux_huntB = B * np.dot(matrix_para(n,par_1)[1],B)
    flux_huntB = flux_huntB[:-1]
    
    flux_I = np.array([flux_I])
    flux_mortB1 = np.array([flux_mortB1])
    flux_huntB1 = np.array([flux_huntB1])
    flux_prelB2 = np.array([flux_prelB2])
    flow = np.concatenate([flux_I,flux_mortB1,flux_huntB1,flux_recy,flux_mortB,flux_prelB2,flux_prelB,flux_huntB])
    return flow




########### EQUILIBRE // EQUILIBRIUM  #################
### Recherche des equilibres // Search for equilibrium  ##########

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver 
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc


# Création du vecteur d'equation // Creation of the equation vector  

def fun (X,n,par_1):
    B1 = X[0]
    B = X[1:]
    dN = par_1["I"] - par_1["l"] * B1 - par_1["a2"]/ par_1["alpha2"] * B1 * B[0] + recyclage(B,n,par_1)[0]
    Y = B * (vect_gi(B1,n,par_1) + np.dot(matrix_para(n,par_1)[0],B) - np.dot(matrix_para(n,par_1)[1],B) - par_1["m"] - par_1["D"] * B)
    dN = np.array([dN])
    W = np.concatenate([dN,Y])
    return W



# Creation de la J estimé // Creation of estimed J (Jacobian)
def jacfct(X,n,par_1):   

    matrix_eq = np.diag(np.ones(n),k=0) 
     
    for i in range (n):
        matrix_eq[i,i] = X[i+1] # matrice des equilibres compartiment biotiques // equilibrium matrice of biotic compartment

    
    A = np.diag(vect_a_b(n,par_1)[1][1:]*vect_a_b(n,par_1)[0][1:]*np.ones(n-1),k = -1) # Matrice d'interaction // Interaction matrix
    A = A + np.diag(-par_1["D"]*np.ones(n),k = 0)
    A = A + np.diag(- vect_a_b(n,par_1)[0][1:]*np.ones(n-1),k = 1)
    
    matrix_fprimeg = vect_gi(X[0],n,par_1) + np.dot(A,X[1:]) - par_1["m"]
    matrix_fprimeg = np.diag(matrix_fprimeg*np.ones(n),k=0) 
    
    
    # Création de la matrice Jacobienne où Jij = Bi %*% dfi(B)/dBj
        
    jac = np.dot(matrix_eq,A) + matrix_fprimeg # exclusivement des compartiments biotiques
    
    jac_add = np.ones(n)*0
    jac_add[0] =  par_1["a2"] * X[1]
    jac_new = np.zeros((n,n+1))
    jac_new[:,0] = jac_add
    jac_new[:,1:] = jac # J_new est la matrice biotique totale comprenant dF(B)/dN
    
    
    jac_new = jac_new.reshape(n*(n+1))
    jac_tt = np.concatenate([J_N(X,n,par_1),jac_new]) # J_N fonction de la matrice jacobienne pour dN/dt
    jac_tt = jac_tt.reshape(n+1,n+1)    # Matrice Jacobienne de mon système totale
    return jac_tt


###### Réponse à la PERTURBATION // Response to perturbation ###############

# Création du la perturbation // Creation of perturbation
    

def sigma_pert(par_1):# Deja la variance de la perturbation chez Shanafelt
    if par_1["sppertu"] == 0:
        sigma_pertu = 0.01
    else:
        sigma_pertu = 0.0025
    
    return sigma_pertu


# Création de la jacobienne partie abiotique DONC ICI B1
# Creation of the Jacobian for abiotic part so B1 // B0 in article

def J_N(X,n,par_1):
    B1 = X[0]
    B = X[1:]
    temp = np.ones(n-1)
    for i in range(n-1):
        temp[i] = (1-vect_a_b(n,par_1)[1][i+1])* vect_a_b(n,par_1)[0][i+1] * vecteur_para(n,par_1)[1][i+1]
    no_conversion_2 = np.diag(temp,k = -1)
    no_conversion_2 = no_conversion_2 + np.diag(temp,k = 1) # use dans le np.dot de Y = matrice Jacobienne du non convertie
    gi_JN = 0 * np.ones(n) # creation du vecteur de croissance/mort independante de la densité
    for i in range (n):
        if i == 0:
            gi_JN[i] = - par_1["a2"]/ par_1["alpha2"]* B1
        else:
            gi_JN[i] = 0
    dN = -par_1["l"] - par_1["a2"]/ par_1["alpha2"]* B[0]
    Y = gi_JN  + par_1["lam"]*vecteur_para(n,par_1)[0] * (par_1["m"] + 2*par_1["D"]*B)  + np.dot(no_conversion_2,B) * par_1["lam"] 
    dN = np.array([dN])
    W = np.concatenate([dN,Y])
    return W


#### Fonction rapport => Savoir si I/Irecymax doit être constant (1 / Oui) ou pas (0 / Non)
# Ratio function => If I/Irecymax must be constant (1 / yes) or not (0 / no)

def rap_i_irecymax(n,par_1):
    if par_1["rapport"] == 1:
        if n == 1:
            if par_1["D"]==0:
                par_1["I_max"] =  0.01
                par_1["I_min"] =  0.01
            else:
                par_1["I_max"] = par_1["l"]/par_1["a2"]**2*(10*par_1["D"]*par_1["l"]+par_1["a2"]*par_1["m"])
                par_1["I_min"] =  par_1["l"]/par_1["a2"]**2*(10*par_1["D"]*par_1["l"]+par_1["a2"]*par_1["m"])
        
        elif n==2:
            if par_1["D"]==0:
                par_1["I_max"] = 0.01
                par_1["I_min"] = 0.01
            else:
                par_1["I_max"] = (par_1["l"]*(10*par_1["a3"]**2*par_1["b3"]*par_1["l"]-par_1["a2"]*par_1["a3"]*par_1["m"]+10*par_1["l"]*par_1["D"]**2+par_1["a2"]*par_1["m"]*par_1["D"]))/(par_1["a2"]**2*par_1["D"])
                par_1["I_min"] = (par_1["l"]*(10*par_1["a3"]**2*par_1["b3"]*par_1["l"]-par_1["a2"]*par_1["a3"]*par_1["m"]+10*par_1["l"]*par_1["D"]**2+par_1["a2"]*par_1["m"]*par_1["D"]))/(par_1["a2"]**2*par_1["D"])
    else:
        par_1["I_max"] =  par_1["I_max"]
        par_1["I_min"] = par_1["I_min"]
    return par_1["I_max"],par_1["I_min"]


