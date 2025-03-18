import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

def Lead_lag_RT(MV,Kp1,Tlag,Tlead,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "FO_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        TRAP: Trapezoïdal method
    
    The function "FO_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp1/(1+K))*((1+Tlead/Ts)*MV[-1]-(Tlead/Ts)*MV[-2]))
            elif method == 'EFD':
                PV.append(((1-K)*PV[-1])+(K*Kp1*((Tlead/Ts)*MV[-1]+(1-(Tlead/Ts))*MV[-2])))
            #elif method == 'TRAP':
                #PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))            
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp1/(1+K))*MV[-1])
    else:
        PV.append(Kp1*MV[-1])
def PID_RT(SP, PV, Man, MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin, MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method='EBD_EBD'):
    #computation of E
    if len(PV) == 0:
        E.append(SP[-1] - PVInit)
    else:
        E.append(SP[-1] - PV[-1])
    #action proportionnelle
    MVP.append(Kc*E[-1])
    #action intégrale
    if len(MVI) == 0:
        MVI.append((Kc*Ts/Ti)*E[-1])
    else :
        if method == "TRAP":
            MVI.append(MVI[-1]+((0.5*Kc*Ts)/Ti)*(E[-1]+E[-2]))
        else : 
            MVI.append(MVI[-1]+((Kc*Ts)/Ti)*E[-1])
    #action dérivée
    Tfd = alpha*Td
    if len(MVD) == 0:
        MVD.append(((Kc*Td)/(Tfd+Ts))*(E[-1]))
    else :
        if method=='TRAP':
            MVD.append((((Tfd-Ts/2)/(Tfd+Ts/2)))*MVD[-1]+((Kc*Td)/(Tfd+Ts/2))*(E[-1]-E[-2]))
        else : 
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+Kc*(Td/(Tfd+Ts))*(E[-1]-E[-2]))
    #integrator reset
    if Man[-1] == True:
        if ManFF:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else:
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]
    #saturation
    if (MVP[-1] + MVI[-1]+MVFF[-1])>MVMax:
        MVI[-1] = MVMax - MVP[-1]-MVFF[-1]
    if (MVP[-1] + MVI[-1]+MVFF[-1])<MVMin:
        MVI[-1] = MVMin - MVP[-1] - MVFF[-1]
   
    MV.append(MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1])
    if MV[-1]>MVMax:
        print(MV[-1])
    
    return (MV, MVP, MVI, MVD)
