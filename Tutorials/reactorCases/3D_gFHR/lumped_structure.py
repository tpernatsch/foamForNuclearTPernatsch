import numpy as np
import matplotlib.pyplot as plt

##compute the different value of heat conductance for the lumped structure.


############### Constant of the problem (radius, thermal conductivities...) #####################
P=280e6 ###power of the enire reactor
N_p=250190 ###number of pebble 
Q_p=P/N_p  ###power of a pebble 
############### Shell constant ###############

R_shell = 2e-2
k_shell = 41.68 #graphite conductivity at 1000K, April thesis

############### Matrix constant ###############

R_m      = 1.8e-2
R_c      = 1.38e-2
k_matrix = 29   # k_graphite=41.68 taking from python effective calculation
rho_m    = 1632 #Novak thesis
Cp_m     = 1698 #Novak thesis at 1000k

############### Trisos constant ###############

R_fuel    = 212.5e-6
R_buffer  = 312.5e-6
R_PyC_in  = 352.5e-6
R_SiC     = 387.5e-6
R_PyC_out = 427.5e-6

k_fuel    = 3.3073
k_buffer  = 0.5 
k_PyC_in  = 4.0
k_SiC     = 90.3 ##Novak thesis at 1000K 
k_PyC_out = 4.0

rho_fuel = 11031
Cp_fuel  = 310 ##Novak thesis at 1000K
V_fuel=4/3*np.pi*R_fuel**3
m_fuel=rho_fuel*V_fuel


rho_SiC=3216
Cp_SiC=1223
V_SiC=4/3*np.pi*(R_SiC**3-R_PyC_in**3)
m_SiC=rho_SiC*V_SiC

rho_PyC=1900
Cp_PyC=720
V_PyC=4/3*np.pi*(R_PyC_out**3-R_buffer**3)-V_SiC
m_PyC=rho_PyC*V_PyC

rho_buffer=1000
Cp_buffer=720
V_buffer=4/3*np.pi*(R_buffer**3-R_fuel**3)
m_buffer =rho_buffer*V_buffer

m_tot=m_fuel+m_buffer+m_SiC+m_PyC
############### Usefull constants ###############

V_p = 4/3*np.pi*R_shell**3 ### pebble volume
V_m=4/3*np.pi*(R_m**3-(R_c)**3)
V_trisos=4/3*np.pi*R_PyC_out**3
N_T = 9022 ### number of trisos iin the matrix

############## Power and density of power ##########

q_matrix=Q_p/V_m ###power density in the matrix
Q_triso=Q_p/N_T  ###Power of one triso 
q_fuel=Q_triso/(4/3*np.pi*R_fuel**3)  ###density of power in the fuel




############### heat conductance between the pebble surface and the average temperature of the matrix ###############

H_shell = (1/(4*np.pi*k_shell)*(1/R_m-1/R_shell))**-1
### unheated matrix
#=(1/(4*np.pi*k_matrix)*(1/(1.64e-2)-1/R_m))**-1 
#if we consider that the matrix is heated
H_sToavM =(1/(V_m*6*k_matrix)*(R_m**2-3/5*(R_m**5-R_c**5)/(R_m**3-R_c**3))+R_c**3/(V_m*3*k_matrix)*(1/R_m-3/2*(R_m**2-R_c**2)/(R_m**3-R_c**3)))**-1
H_3=H_sToavM*H_shell/(H_sToavM+H_shell)

############### heat conductance between the average matrix temperature and the average fuel temperature 

H_avMTosf = (1/(4*np.pi*N_T)*(1/k_buffer*(1/R_fuel-1/R_buffer)+1/k_PyC_in*(1/R_buffer-1/R_PyC_in)+1/k_SiC*(1/R_PyC_in-1/R_SiC)+1/k_PyC_out*(1/R_SiC-1/R_PyC_out)))**-1
H_sfToavf= ((1/(4*np.pi*R_fuel**3*N_T))*(R_fuel**2/(15*k_fuel)))**-1
H_2=H_avMTosf*H_sfToavf/(H_avMTosf+H_sfToavf)

############### heat conductance between the average temperature of the fuel and the maximum temperature of the fuel

H_sfTomaxf = (1/(4*np.pi*R_fuel**3*N_T)*R_fuel**2/(6*k_fuel))**-1
H_1        = 10/15*H_sfToavf
#H_sfTomaxf*H_sfToavf/(H_sfToavf-H_sfTomaxf)


############### first node at the avreage matrix temperature ###############
rho_2 = rho_m
Cp_2  = Cp_m

############### second node at the average fuel temperature ###############

Cp_1=(Cp_fuel*m_fuel+Cp_buffer*m_buffer+Cp_SiC*m_SiC+Cp_PyC*m_PyC)/m_tot
rho_1=m_tot/V_trisos

############### Volume of each node ###############
V_1 = V_trisos
### Volume of the matrix
V_3=(4/3*np.pi*R_shell**3-4/3*np.pi*R_m**3)
### Volume of the pebble
V_tot=V_p
############### Volumetric area ###############
V_fuel=4/3*np.pi*R_fuel**3
V_area=3*0.61/2e-2

############### Power fraction ######################


###take into account the fact only pebbles are producing power in the reactor 
q_real=Q_p/V_p 

##take into account the volume fraction lower than 1.

V_frac_real=V_trisos*N_T/V_tot+(V_3+V_m-V_trisos*N_T)/V_tot

############### Print results ###############


print('heat conductance 1, 2 and 3, from the maximum to the surface',H_1/V_p,H_2/V_p,H_3/V_p,sep='  ')
print('rhoCp for each node, from the max to the surface', rho_1*Cp_1, rho_2*Cp_2/V_frac_real, sep='  ')
print('fraction of volume for each node, frome the max to the surface',V_trisos*N_T/V_tot,(V_3+V_m-N_T*V_trisos)/V_tot,sep='  ')
print('real power density',q_real,sep=':')


