# -*- coding: utf-8 -*-
"""
# 
Devoir 2
Created on Mon Oct 10 09:53:01 2022
Ecole Polytechnique Montreal - Universite de Montreal
MEC 8211 : Vérif. et valid. en modélisation numérique
@auteurs : Houssem Eddine Younes & Mohamed Dhia Slama
October 2022

"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.sparse.linalg.dsolve import linsolve 


def terme_source(T_MMS,D_eff,k):
     x,y = sp.symbols('x y')
     
     #Application de l'opérateur sur la solution MMS
     source = sp.diff(T_MMS,y)-D_eff*(1/x)*(sp.diff(T_MMS,x)+x*sp.diff(sp.diff(T_MMS,x),x))+k*T_MMS
     
     #Création d'une fonction symbolique
     f_source_MMS = sp.lambdify([x,y], source, "numpy")
     
     
     #Derivee partielle pour neumann 
     deriv_neumann = sp.diff(T_MMS,x)
     f_neumann = sp.lambdify([x,y], deriv_neumann, "numpy")
     
     
     return f_source_MMS,f_neumann
def calculations_kernel(i_schemas): 
    #Initialisation
    Matrix = np.zeros([N_r+1,N_r+1])
    C_new = np.zeros([N_r+1])
    C_old = np.zeros([N_r+1])
    r=np.linspace(0,R,N_r+1)
    
    #Terme source 
    f_source_MMS,f_neumann = terme_source(T_MMS,D_eff,k)   
    
    #Creation de la matrice A 
    #Discretisation de gear pour la CL de neumann a r = 0 
    Matrix[0,0] = -3
    Matrix[0,1] = 4 
    Matrix[0,2] = -1 
    Matrix[N_r,N_r] = 1
    
    for i in range(1,N_r): 
        A,B,C = schemas_numeriques(i,i_schemas,r)
        Matrix[i,i] =  A 
        Matrix[i,i-1] = C 
        Matrix[i,i+1] = B 
    
    #Conditions initiales 
    sol_MMS = sp.lambdify([x,y], T_MMS, "numpy")
    C_old[:] = sol_MMS(r[:],0) 
    
    
    #Ajout du terme source dans le vecteur de droite 
    C_old[1:N_r] -= f_source_MMS(r[1:N_r],0)
    
    #Conditions aux limites pour la premiere iteration 
    C_old = boundary_conditions(Matrix,C_old,0,T_MMS,0,f_neumann)
    # C_old[N_r] = Ce
    
    
    C_test = np.copy(C_old)
    
    #initialisation
    Erreur = 100
    t= 0 
    DT_r=0
    iterations = 0 
    Boolean = False
    #Boucle iterative avec une condition de convergence temporel
    while Boolean == False:
        iterations += 1 
        t = t+1
        DT_r = DT_r+ Delt_t
        C_new = np.linalg.solve(Matrix,C_old)
        #Calcul d'erreur entre iteration : Condition d'arret
        Erreur = sum(abs(C_new-C_test))
        #Assignation de C caluclee a C_old
        C_old = np.copy(C_new)
        C_test = np.copy(C_old)
        #Terme source
           
        C_old[1:N_r] = C_old[1:N_r] -f_source_MMS(r[1:N_r],DT_r)
        f_source_vect = f_source_MMS(r[1:N_r],DT_r)
            
        C_old = boundary_conditions(Matrix,C_old,0,T_MMS,DT_r,f_neumann)
        if Erreur <= 10**(-4) : 
            Boolean = True
        elif iterations > 50 : 
            Boolean = True
        print("iteration "+str(t)+"," + "erreur :" + str(Erreur))
        
    
    kkkk = 0 
    return C_new,r,DT_r,f_neumann
def boundary_conditions(Matrix,C_old,index,T_MMS,DT_r,f_neumann): 
    sol_MMS = sp.lambdify([x,y], T_MMS, "numpy")
        
    #Dirichlet pour r = R
    C_old[N_r] = sol_MMS(R,DT_r)
    #Neumann
    C_old[0] = f_neumann(0,DT_r)*2*Delt_r
    ll = f_neumann(0,DT_r)*2*Delt_r
    
    return C_old       
def schemas_numeriques(i,schema,r): 
    # Schemas avec termes source constant (pour le devoir 1)
    if schema == 0 :  # 1er schema
        #Diag
        A = (1+D_eff*Delt_t/(Delt_r*r[i])+2*D_eff*Delt_t/(Delt_r**2))
        #Ci+1
        B =  (-D_eff*Delt_t/(Delt_r*r[i])-D_eff*Delt_t/(Delt_r**2))
        #Ci-1
        C= (-D_eff*Delt_t)/(Delt_r**2)
    elif schema == 1: # 2eme schema
        # #Diag
        A = (1+2*D_eff*Delt_t/(Delt_r*Delt_r))
        #Ci+1
        B = -(D_eff*Delt_t/(Delt_r*Delt_r))-(D_eff*Delt_t/(2*r[i]*Delt_r))
        #Ci-1
        C = (D_eff*Delt_t/(2*r[i]*Delt_r))-(D_eff*Delt_t/(Delt_r*Delt_r))
        
    #Schemas avec terme source de premier ordre (pour le devoir 2 ) question D
    elif schema ==2 : # 1er schema VAR
        #Diag
        # A = (Delt_r**2+D_eff*Delt_t*Delt_r/r[i]+2*D_eff*Delt_t)/((Delt_r**2)*(1-k*Delt_t))
        A = (Delt_r*Delt_r+D_eff*Delt_t*Delt_r/r[i]+2*D_eff*Delt_t)/(Delt_r*Delt_r*(1-k*Delt_t))
        #Ci+1
        C=(-D_eff*Delt_t)/(Delt_r*Delt_r*(1-k*Delt_t))
        # B= (-D_eff*Delt_t*Delt_r/r[i]-D_eff*Delt_t)/((Delt_r**2)*(1-k*Delt_t))
        #Ci-1
        B = (-D_eff*Delt_t*Delt_r/r[i]-D_eff*Delt_t)/(Delt_r*Delt_r*(1-k*Delt_t))
        # C= (-D_eff*Delt_t)/((Delt_r**2)*(1-k*Delt_t))
    elif schema ==3: #2eme schema VAR
        #Diag
        A = (Delt_r**2+2*D_eff*Delt_t) /((Delt_r**2)*(1-k*Delt_t))
        #Ci+1
        # B =  (-D_eff*Delt_r*Delt_t/(2*r[i]-D_eff*Delt_t))/((Delt_r**2)*(1-k*Delt_t))
        C= (D_eff*Delt_t*Delt_r/(2*r[i])-D_eff*Delt_t)/(Delt_r*Delt_r*(1-k*Delt_t))
        #Ci-1
        # C = (D_eff*Delt_t*Delt_r/(2*r[i])-D_eff*Delt_t)/((Delt_r**2)*(1-k*Delt_t))
        B = (-D_eff*Delt_t*Delt_r/(2*r[i])-D_eff*Delt_t)/(Delt_r*Delt_r*(1-k*Delt_t))
    return A,B,C
    


def Plot_anal_num(C_new,iterations,r,nb_iter,T_MMS,DT_r,f_neumann): 
   
    
    x, y = sp.symbols('x y')
    # T_analy= 0.25*(S/D_eff)*(R**2)*(x*x/(R**2)-1)+Ce
    # sol_analytique_int = sp.lambdify([x], T_analy, "numpy")
    # sol_analytique = sol_analytique_int(r)
    sol_analytique_int = sp.lambdify([x,y], T_MMS, "numpy")
    temps_fin=np.linspace(0,R,N_r+1)
    temps_fin[:] = DT_r #10**2 # a revoir peut etre il faut mettre DT_r
    #ri, ti = np.meshgrid(r,temps_fin)
    
    sol_analytique = sol_analytique_int(r,temps_fin)
    #sol_analytique[0] = f_neumann(0,DT_r)*2*Delt_r
    plt.plot(r,sol_analytique, label="Solution analytique")
    
    
    plt.plot(r,C_new,'bo' ,label="Solution numerique")
    plt.xlim(0,R)
    y_titre="Concentration molaire"
   
    #plt.title(label="Coupe en y=" + y_titre +  " pour " + str(number_of_elements) + " éléments (" + type_element[self.Type_maillage_index] +")" + " avec le " +type_solveur[schema]  + " (Peclet =" + str(self.Peclet[self.Index_Peclet]) + ")" ,fontsize=11,color="black")
    plt.xlabel("Rayon r [m]")
    plt.ylabel("C , Concentration en sel [mol/m3]")
    String = "Solution analytique vs solution numerique avec "  + str(N_r) + " noeuds" + " pour le schema " + str(i_schemas+1)
    plt.title(label=String ,fontsize=11,color="black")
    plt.legend()
    plt.show()
    
    return sol_analytique


def calcul_erreur(nb_iter,C_new,N_t,sol_analytique): 
    
    L1error[nb_iter] = (1/N_r) * np.sum(np.abs(C_new-sol_analytique))
    # L2error[nb_iter] =  np.sqrt((1/R)*np.sum(N_r*np.square(C_new-sol_analytique)))
    L2error[nb_iter] =  np.sqrt((1/N_r)*np.sum(np.square(C_new-sol_analytique)))
    # L2error[nb_iter] =  np.sqrt(1/N_r*Delt_r*np.sum(np.square(C_new-sol_analytique)))
    Linf_error[nb_iter] = np.max(np.abs(C_new-sol_analytique))
                                  
def plot_errors (): 
    #Plotting Ln(E) vs Ln(1/nx)
    
    # error = 0.5*iterations**2
    # a,b = np.polyfit(np.log(iterations), np.log(L2error), 1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    #errors plotting
    ax.plot(1/iterations,L2error,'-o',label ="L2error")
    ax.plot(1/iterations,L1error,'-o',label ="L1error")
    ax.plot(1/iterations,Linf_error,'-ko',label ="Linf_error")
   
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(b=True, which='minor', color='grey', linestyle='--')
    ax.grid(b=True, which='major', color='k', linestyle='--')
    plt.ylabel('$\Vert erreur \Vert $')
    plt.xlabel('$\ 1/Nr $')
    plt.title(label='Courbes log-log des erreurs pour le schema ' +str(i_schemas+1) ,fontsize=14,color="black")
    plt.legend(loc="upper left")
    plt.figtext(.01, .02, "Ordre P = " +str(np.abs(np.round(Order[len(Order)-1],5))))
    plt.show()
        
##########################################################################################################################
#                                                                                                                       ##
#                                                                                                                       ##
#                                              BOUCLE PRINCIPALE                                                        ##
#                                                                                                                       ##
##########################################################################################################################
#Input Data 
D = 1 # en m 
R = D/2

#Coefficient de diffusion effectif 
D_eff = 10**(-10) # en m2/s

#Constante de reaction 
k = 2 * 10**(-9) # en s-1

#Concentration de sel 
Ce = 10 # mol/m3

#Terme source 
S = 10**(-8) #mol/m3/s

#Solution MMS
#Solution analytique
c_zero = 3.
t_zero1 = 10.
x,y = sp.symbols('x y')
T_MMS = c_zero*sp.cos(y/t_zero1)*sp.sin(np.pi*x/R)

#Definition des variables pour les erreurs et ordre 
# iterations = np.array([20,40,80,160,320])
# iterations = np.array([160,320,640,1280,2560])
iterations = np.array([5,10,20,40,80])
#iterations = np.array([3,6,12,24,48])
L1error = np.zeros(len(iterations-1))
L2error = np.zeros(len(iterations-1))
Linf_error = np.zeros(len(iterations-1))
Order= np.zeros(len(iterations)-1)

#Boucle sur les differents schemas 
for i_schemas in range(3,4): 
    
    #Boucle de raffinement
    for nb_iter in range(len(iterations)):
    
        #Variables de l'espace 
        N_r = iterations[nb_iter] 
        Delt_r = (R)/N_r

        #Variables de temps
        T_max = 100 # en s 
        N_t = 10000
        Delt_t = T_max/N_t
    
        C_new,r,DT_r,f_neumann = calculations_kernel(i_schemas)
        sol_analyt = Plot_anal_num(C_new,iterations,r,nb_iter,T_MMS,DT_r,f_neumann)
        calcul_erreur(nb_iter,C_new,N_t,sol_analyt) 
    
        if nb_iter > 0 : 
            #Calcul d'ordre seulement pour la norme L2
            Order[nb_iter-1] = np.log(L2error[nb_iter]/L2error[nb_iter-1])/np.log(2)
      
    plot_errors()    