# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:45:01 2020

@author: Theis Kevin
"""
from fonction_3_0 import* 
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.optimize import root
import random as rdn
import multiprocessing as mp

def main(line):
    
    par_1 = dictio[line]
    n = int(par_1["n"])
    
    # Ratio function => in parameters if "rapport" = 0, i is constant, if "rapport" = 1, i varies as a function of D
    # In article rapport = 0 for all 
    rap_i_irecymax(n,par_1) # i cst si rapport = 0 / i varie en fct de D (i/irecymax = cst) si r rapport = 1


    
    nb_iteration = int((dictio[line]["lam_max"]-dictio[line]["lam_min"])/dictio[line]["lam_step"] + 1)
    lam_max = dictio[line]["lam_max"]
    lam_min = dictio[line]["lam_min"]
    lam = np.linspace(lam_max,lam_min, int(nb_iteration) )  # nb_iteration
    
    nb_I = int((dictio[line]["I_max"]-dictio[line]["I_min"])/dictio[line]["I_step"] + 1)
    I_max = dictio[line]["I_max"]
    I_min = dictio[line]["I_min"]
    seq_I = np.linspace(I_min,I_max, int(nb_I) )  # nb_iteration
    
    # Matrice qui va recevoir toutes les données // Result Matrix that get all data
    matrix_result = []
    
    
    for inp in seq_I:
        par_1["I"] = inp
        par_1["I_old"] = par_1["I"] # mise en memoire de la valeur d'input d'origine
        
        for lamda in lam:
            par_1["lam"] = lamda
            print(par_1["lam"])
            print("valeur de I avt compen : ", par_1["I"])
           
            
            
            ## Valeur de Irecymax pour TL1 car on commence a lambda = 0.9
            if n == 1:          
                if lamda == max(lam): 
                    par_Irecy = dict(par_1)
                    par_Irecy["lam"] = 1
                    equi_lam = np.ones(n)
                    if par_Irecy["D"] == 0:
                        equi_lam[0] = (10 * par_Irecy["I"])/ par_Irecy["m"]
                    else:
                        equi_lam[0] =(par_Irecy["I"]*par_Irecy["a2"]-par_Irecy["m"]*par_Irecy["l"])/(par_Irecy["D"]*par_Irecy["l"])
                I_maxrecy = recyclage(equi_lam,n,par_Irecy)[0]
#                print("Irecymax1 =",I_maxrecy)
            
                
        ## ZONE de compensation de la perte 
            if par_1["compen"] == 1:
                if lamda != 1: # afin d'avoir un input constant entre le recyclage et non recyclage
                    par_1["I"] = par_1["I_old"]
                    par_1["I"] = par_1["I"] + (1-lamda) * I_maxrecy
#                    print("lam!=", lamda)
                    print("new input", par_1["I"])
                
                

            
            # Methode par integration pour trouver les equilibres // Integration part in order to find equilibrium
            
            tfinal = 4000
            dt = 1e-2 #pas de temps 
            Tps = np.linspace(0, tfinal, int(tfinal/dt+1))
            condi_ini = np.ones(n+1)
            if n ==1:
                if par_1["D"] == 0:
                    for condi in range(n+1):
                        if condi == 0:
                            condi_ini[condi] = par_1["m"]/par_1["a2"]
                        else:
                            condi_ini[condi] = 2
                    
                else:
                    for condi in range(n+1):
                        if condi == 0:
                            condi_ini[condi] = par_1["I"]/par_1["l"] 
                        else:
                            condi_ini[condi] = 10
            elif n ==2:
                if par_1["D"] == 0:
                    for condi in range(n+1):
                        if condi == 0:
                            condi_ini[condi] = par_1["I"]/par_1["l"] 
                        elif condi ==1:
                            condi_ini[condi] = par_1["m"]/(par_1["a3"]*par_1["b3"]) 
                        else:
                            condi_ini[condi] = (par_1["a2"]*par_1["I"]-par_1["l"]*par_1["m"])/(par_1["a3"]*par_1["l"])+1e-5
                    
                else:
                    for condi in range(n+1):
                        if condi == 0:
                            condi_ini[condi] = par_1["I"]/par_1["l"] + 5
                        elif condi ==1:
                            condi_ini[condi] = (par_1["a2"]*par_1["D"]*par_1["I"]+par_1["a3"]*par_1["l"]*par_1["m"]-par_1["D"]*par_1["l"]*par_1["m"])/(par_1["a3"]**2*par_1["b3"]*par_1["l"]+par_1["D"]**2*par_1["l"]) 
                        else:
                            condi_ini[condi] = (par_1["a2"]*par_1["b3"]*par_1["a3"]*par_1["I"]-par_1["a3"]*par_1["b3"]*par_1["l"]*par_1["m"]-par_1["D"]*par_1["l"]*par_1["m"])/(par_1["a3"]**2*par_1["b3"]*par_1["l"]+par_1["D"]**2*par_1["l"])+1e-5
                
            else:
                for condi in range(n+1):
                        if condi == 0:
                            condi_ini[condi] = par_1["I"]/par_1["l"] +1
                        elif condi ==1:
                            condi_ini[condi] = (par_1["a2"]*par_1["I"]*(par_1["a4"]**2*par_1["b4"]+par_1["D"]**2)- par_1["l"]*par_1["m"]*(par_1["a4"]**2*par_1["b4"]+par_1["a3"]*par_1["a4"]-par_1["D"]*par_1["a3"]+par_1["D"]**2))/(par_1["l"]*par_1["D"]*(par_1["a3"]**2*par_1["b3"]+par_1["a4"]**2*par_1["b4"]+par_1["D"]**2)) +1
                        elif condi == 2:
                            condi_ini[condi] = (par_1["a2"]*par_1["a3"]*par_1["b3"]*par_1["I"]+par_1["l"]*par_1["m"]*(par_1["a4"]-par_1["a3"]*par_1["b3"]-par_1["D"]))/(par_1["l"]*(par_1["a3"]**2*par_1["b3"]+par_1["a4"]**2*par_1["b4"]+par_1["D"]**2)) + 1
                        else:
                            condi_ini[condi] = (par_1["a2"]*par_1["a3"]*par_1["a4"]*par_1["b3"]*par_1["I"]*par_1["b4"]-par_1["l"]*par_1["m"]*(par_1["a3"]**2*par_1["b3"]+par_1["a3"]*par_1["a4"]*par_1["b3"]*par_1["b4"]+par_1["D"]*par_1["a4"]*par_1["b4"]+par_1["D"]**2))/(par_1["l"]*par_1["D"]*(par_1["a3"]**2*par_1["b3"]+par_1["a4"]**2*par_1["b4"]+par_1["D"]**2))+1
                            
                            
            
            X = odeint(syst,condi_ini,Tps,args=(n, par_1))
            Eq = X[-1]
    
            ## Permet d'ajuster l'intégration si besoin d'un temps d'inte plus long
            # Allows you to adjust the integration if you need more time to integrate
            
            if np.all(np.abs(Eq - X[-5]) <= 10**-5) == False : # Si valeur négative, mets un message d'erreur !
                tfinal = 40000
                dt = 1e-2 #pas de temps
                Tps = np.linspace(0, tfinal, int(tfinal/dt+1))
                condi_ini = condi_ini
            
                X = odeint(syst,condi_ini,Tps,args=(n, par_1))
                Eq = X[-1]
                print("changement d'EQ 1")
                
                if np.all(np.abs(Eq - X[-5]) <= 10**-5) == False : # Si valeur négative, mets un message d'erreur !
                    tfinal = 1000000
                    dt = 1e-2 #pas de temps
                    Tps = np.linspace(0, tfinal, int(tfinal/dt+1))
                    condi_ini = condi_ini
                
                    X = odeint(syst,condi_ini,Tps,args=(n, par_1))
                    Eq = X[-1]
                    print("changement d'EQ 2 ")     
                    

                    if np.all(np.abs(Eq - X[-5]) <= 10**-5) == False :
                        raise Exception ("Pas atteint le point d'équilibre")
                    
            if np.any(Eq<0) == True :  # si un equilibre est inférieur a zéro remet a zero
                for i in range(len(Eq)):
                    if Eq[i]<0:
                        Eq[i] = 0
                        
            print ("Les équilibres sont", Eq)
                    
                    
                    
                    
            
    
            ###### Réponse à la PERTURBATION // Response to PERTURBATION ###############
    
            
           ########Trouver la matrice de COVARIANCE et de CORRELATION // Find the COVARIANCE and CORRELATION matrix  ########################
            
            #### Façon 1 via les produits de Kronecker mais pas encore au point 
            # Va définir le type de perturbation
            
            T = np.diag(np.ones(n+1),k=0) 
         
            for i in range (n+1):
                if i == 0:
                    T[i,i] = Eq[i]**(par_1["typepertuB1"]) # Perturbation des nutriments types environnementales si != ** 1
                else :
                    T[i,i] = Eq[i]**(par_1["typepertuB"])
                    # T le type de perturbation mise en place sur les éléments ici demographique
            T_add = np.zeros(n+1)
            T_add[0] = par_1["I"]
            T_new = np.zeros((n+1,n+2))
            T_new[:,0] = T_add
            T_new[:,1:] = T # Matrice des perturbations en considerant I / N / Bi = as Shanafelt article
            
            
            # Matrice variance-covariance de la stochasticité
            Ve = np.diag(np.ones(n+2)*0,k = 0)
            # NB : n+1 car on ajoute I (afflux de N) et N en plus de B
            par_1["sppertu"] = int(par_1["sppertu"])
            Ve[par_1["sppertu"],par_1["sppertu"]] = sigma_pert(par_1)
            
            
            ##### PARTIE creation de la J // Part : Creation of Jacobian #####
            
            matrix_eq = np.diag(np.ones(n),k=0) 
    
            for i in range (n):
                matrix_eq[i,i] = Eq[i+1] # matrice des equilibres compartiment biotiques // Equilibrium matrix of biotic compartment
    
            
            A = np.diag(vect_a_b(n,par_1)[1][1:]*vect_a_b(n,par_1)[0][1:]*np.ones(n-1),k = -1) # Matrice d'interaction // Interaction matrix 
            A = A + np.diag(-par_1["D"]*np.ones(n),k = 0)
            A = A + np.diag(- vect_a_b(n,par_1)[0][1:]*np.ones(n-1),k = 1)
            
            matrix_fprimeg = vect_gi(Eq[0],n,par_1) + np.dot(A,Eq[1:])- par_1["m"]
            matrix_fprimeg = np.diag(matrix_fprimeg*np.ones(n),k=0)
            
            # Création de la matrice Jacobienne où Jij = Bi %*% dfi(B)/dBj
            
            # (fg)' = f'g + fg'
            
            J = np.dot(matrix_eq,A) + matrix_fprimeg # exclusivement des compartiments biotiques // partie => fg' + f'g
            
            ### Calcule a part des d²Bi/dtdBi afin d'ajouter le gamma qui represente la force de le self regulation
            
            
            
            J_add = np.ones(n)*0
            J_add[0] =  par_1["a2"] * Eq[1]
            J_new = np.zeros((n,n+1))
            J_new[:,0] = J_add
            J_new[:,1:] = J # J_new est la matrice biotique totale comprenant dF(B)/dN
            
            
            J_new = J_new.reshape(n*(n+1))
            J_tt = np.concatenate([J_N(Eq,n,par_1),J_new]) # J_N fonction de la matrice jacobienne pour dN/dt
            matrice_jacobian = J_tt
            J_tt = J_tt.reshape(n+1,n+1)    # Matrice Jacobienne de mon système totale // Jacobian matix of all system
            
            
            J_inverse = np.linalg.inv(J_tt)
            
            J_test = np.dot(J_inverse,J_tt)
            
            J_inverse = J_inverse.reshape((n+1)*(n+1)) 
            
            
            matrix_compa = np.diag(np.ones(n+1)*1,k = 0)
            
            if np.allclose(J_test,matrix_compa) == False:
                raise Exception ("Jacobian inverse pas bonne")
           
            
            
        # Valeur propre dominante // Dominant eigenvalue
        
            VP = np.linalg.eigvals(J_tt)
    
            VP = max(np.real(VP))
    
            
            # Resolution afin de trouver la matrice de variance-covariance des variables (des interactios entre sp)
            
            I_matrix = np.diag(np.ones(n+1),k = 0)
            temp1 = np.kron(J_tt,I_matrix) + np.kron(I_matrix,J_tt)
            temp = - np.linalg.inv(temp1)
            TVT = np.dot(T_new,Ve)
            TVT = np.dot(TVT,np.transpose(T_new))
            TVT = np.reshape(np.transpose(TVT),(n+1)**2)
            C_star_line = np.dot(temp,TVT)
            C_star_matrix = C_star_line.reshape(n+1,n+1) # matrice de variance-covariance des interactios entre sp
            
            
            R_star_matrix = np.diag(np.ones(n+1)*4,k = 0) # matrice de correlation // Correlation matrix
            for i in range(n+1):
                for j in range (n+1):
                    R_star_matrix[i,j] = C_star_matrix[i,j]/(C_star_matrix[i,i]*C_star_matrix[j,j])**0.5
            
            R_star_line = R_star_matrix.reshape((n+1)*(n+1))
            flux_Eq = flux(Eq,1,n,par_1)


            # Ajout des valeurs de recyclages avec les equilibres Eq
            
            entre_recy = Eq[1:]
            valu_recy = recyclage(entre_recy,n,par_1)
            print("recy tot:", valu_recy[0])
            
            # ZONE de compensation de la perte // Loss Compensation Area
            if n != 1 :
                if lamda == max(lam): # petite triche ici attention 
                    par_Irecy = dict(par_1)
                    par_Irecy["lam"] = 1
                    equi_lam = np.ones(n) # entre_recy
                    if n == 2:
                        if par_Irecy["D"] == 0:
                            equi_lam[0] = entre_recy[0]
                            equi_lam[1] = entre_recy[1]
                        else:
                            equi_lam[0] = entre_recy[0]
                            equi_lam[1] = entre_recy[1]
                    if n == 3:
                        if par_Irecy["D"] == 0:
                            equi_lam[0] = entre_recy[0]
                            equi_lam[1] = entre_recy[1]
                            equi_lam[2] = entre_recy[2]
                        else:
                            equi_lam[0] = entre_recy[0]
                            equi_lam[1] = entre_recy[1]                
                            equi_lam[2] = entre_recy[2]
                            
                    I_maxrecy = recyclage(equi_lam,n,par_Irecy)[0]

                
            par_1["I_Irecymax"] = par_1["I"]/I_maxrecy



            # Traitement de resultat // Result treatment
            
            # RAJOUT par_1["I_Irecymax"] entre I et l
            
            val_para = [par_1["rapport"], par_1["compen"], par_1["n"] , par_1["typepertuB1"] , par_1["typepertuB"] , par_1["sppertu"], par_1["I_old"] , par_1["I"], par_1["I_Irecymax"] , par_1["l"] , lamda , par_1["D"], par_1["a2"], par_1["a3"], par_1["a4"], par_1["a5"], par_1["b2"], par_1["b3"], par_1["b4"], par_1["b5"],par_1["alpha3"], par_1["m"] ]
            VP = [VP]
            
            
            pre_result = np.concatenate([val_para,Eq,valu_recy,C_star_line,R_star_line,VP , flux_Eq, matrice_jacobian, J_inverse])
            
            
            matrix_result.append(pre_result)
            
            
    return matrix_result


if __name__ == '__main__' :
    print("avt",len(dictio))
    i = range(len(dictio))
    p = mp.Pool(4)
    result = p.map(main,i)
    valu_n = []
    for k in i:
        valu_n.append(dictio[k]["n"])
     
    n = int(max(valu_n))
    print("après",n)
    
    
    # Création du ficher // File creation
    
    filename = "data_D_TL1_pertuB0.txt"
    file = open(filename,"w+")
    
    #####  /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
    ### ATTENTION si lambda différent il faut regarder si par_recymax correspond bien !!!!!
    
    string_para ="rapport" + ";"+ "compen" + ";" + "n" + ";" + "type_pertuB1" + ";" + "type_pertuB" + ";" + "sp_pertu" + ";" + "I_old" + ";" + "I" + ";" + "I_Irecymax" + ";" + "lessi" + ";" + "lambda" +  ";" + "D" + ";" + "a2" + ";" + "a3" + ";" + "a4" + ";" + "a5" + ";" + "b2" + ";" + "b3" + ";" + "b4" + ";" + "b5" + ";" + "alpha3" + ";" + "m" 
    
    ls_eq = [] # liste des headers a d'equilibre
    ls_recy = ["recy_tot","recy_mort","recy_no_convertie",] # liste des headers pour le recyclage
    
    
    for l in range(n+1):
        ls_eq.append(("Eq_B"+str(l+1)))

    
    
    for l in range(1,n+1):
        ls_recy.append(("recy_B"+str(l+1)))
        
    for l in range(1,n+1):
        ls_recy.append(("recy_mort_B"+str(l+1)))   
    for l in range(1,n+1):
        ls_recy.append(("recy_noconvers_B"+str(l+1)))      
        
        
    
    string_eq = str()
    string_recy = str()
    
    for l in range (len(ls_eq)):
        string_eq = string_eq + ";" + ls_eq[l]
    
    for l in range (len(ls_recy)):
        string_recy = string_recy + ";" + ls_recy[l]    
    
    ls_VC = [] # liste des headers de VC (matrice de variance/covariance)
    ls_R = [] # liste des headers de R (matrice de correlation)
    
    for l1 in range (n+1):
        for l2 in range(n+1):
            ls_VC.append(("VC_"+str(l2+1)+str(l1+1)))
            ls_R.append(("R_"+str(l2+1)+str(l1+1)))

            
    
    string_VC = str()
    string_R = str()

    
    for l in range (len(ls_VC)):
        string_VC = string_VC + ";" + ls_VC[l]
        string_R = string_R + ";" + ls_R[l]
   
        
    string_VP = ";" + "VP"
    

    ls_flux_Eq = ["fluxEq_I","fluxEq_mortB1","fluxEq_huntB1",]
    
    for l in range(1,n+1):
        ls_flux_Eq.append(("fluxEq_recyB"+str(l+1)))
        
    for l in range(1,n+1) :  
        ls_flux_Eq.append(("fluxEq_mortB"+str(l+1))) 
        
    for l in range(1,n+1):
        ls_flux_Eq.append(("fluxEq_prelB"+str(l+1))) 
        
    for l in range(1,n):
        ls_flux_Eq.append(("fluxEq_huntB"+str(l+1))) 
        

    string_fluxEq = str()
    
    for l in range (len(ls_flux_Eq)):
        string_fluxEq = string_fluxEq + ";" + ls_flux_Eq[l]
        
    

    ls_J = []
    ls_Jinverse = []    
        
    for l in range(n+1):
        for m in range(n+1):
            ls_J.append("J"+str(l+1)+str(m+1))
            ls_Jinverse.append("J_inverse"+str(l+1)+str(m+1))
        
    string_J = str()
    string_Jinverse = str()
    
    for l in range (len(ls_J)):
        string_J = string_J + ";" + ls_J[l]
        string_Jinverse = string_Jinverse + ";" + ls_Jinverse[l]
        
    
        
    file.write(string_para + string_eq + string_recy  +string_VC +  string_R + string_VP + string_fluxEq + string_J + string_Jinverse + ";" + "supp" + "\n")
    
    
     # Donnée en train d'être écrit
    
    strs = [ "" for x in range(len(result[0])*len(result))] # creation de mon nombre de string (1 par type de paramètre)
    


    for i in range(len(result)):
        for j in range(len(result[0])):
            for niv_3 in range(len(result[0][0])):
                strs[i*len(result[0])+j] = strs[i*len(result[0])+j] + str(result[i][j][niv_3])+ ";"
    
    
    for i in range(len(strs)):
        file.write(strs[i] +  "\n")

    file.close()         