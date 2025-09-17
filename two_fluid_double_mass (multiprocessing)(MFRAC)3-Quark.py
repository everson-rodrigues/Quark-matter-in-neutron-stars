#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 14:58:31 2024

@author: everson
"""

import multiprocessing
import numpy as np
import pandas as pnd
import scipy.integrate as spi
from scipy import interpolate
import matplotlib.pyplot as plt
import time
import csv
import gc

start_time = time.time()
# Constants
gms = 1.475562  # [L]/[M]
ms = 5660.57
h = 197.3269631

nme_frm = []
dme_frm = []
rnm_frm = []
rdm_frm = []
totalmass_frm = []
totalradius_frm = []
starradius_frm = []
list_fr = []
index_fr = []
deltaradius_frm = []


nstars = 8000
nsteps = 1
niterations = 20
search_fr=1e-5

frm = 0.12
# fr=0.8
fr=0.005
frf=0.105

prof_mass1=0.8
prof_mass2=1.0
prof_mass3=1.4
prof_mass4=1.9

while fr<=frf:
    frstr=str("{:.3f}".format(fr))
    print("---------------------------------------------------------------------")

    print(frstr)
    file_name=f'youreos.txt'
    print(file_name)
    
    data = pnd.read_csv(file_name,dtype=None, delim_whitespace=True,header=None)
    # data = data[::-1]
    rho=data.iloc[:,0]
    
    eNM=data.iloc[:,1]
    pNM=data.iloc[:,2]
    
    eDM=data.iloc[:,3]
    pDM=abs(data.iloc[:,4])
    
    
    
    # Interpolate the pressure as function of the energy densitye and vice-versa
    efpNM=interpolate.interp1d(pNM,eNM,fill_value='extrapolate', kind='linear')
    efpDM=interpolate.interp1d(pDM,eDM,fill_value='extrapolate', kind='linear')
    
    pfeNM=interpolate.interp1d(eNM,pNM,fill_value='extrapolate', kind='linear')
    pfeDM=interpolate.interp1d(eDM,pDM,fill_value='extrapolate', kind='linear')
    
    # Data Stars 
    eNM_datastars=np.exp(np.linspace(np.log(eNM[1]), np.log(eNM[eNM.shape[0]-1]), nstars))
    eDM_datastars=np.exp(np.linspace(np.log(eDM[1]), np.log(eDM[eDM.shape[0]-1]), nstars))
    # eNM_datastars=np.linspace(eNM[1], eNM[eNM.shape[0]-1], nstars)
    # eDM_datastars=np.linspace(eDM[1], eDM[eDM.shape[0]-1], nstars)
    
    pNM_datastars=pfeNM(eNM_datastars)
    pDM_datastars=pfeDM(eDM_datastars)
    
    
    # TOV equations for two fluids
    def tov_two_fluids(r,y):
        pn, pd, mn, md = y  # m = total mass, pn = pressure for normal matter, pd = pressure for dark matter
            
        dpndr = - (pn+efpNM(pn)) * ((mn+md) + 4 * np.pi * r**3 *(pn+pd)/ms)/(r**2/gms-2*r*(mn+md))
        dpddr = - (pd+efpDM(pd)) * ((mn+md) + 4 * np.pi * r**3 *(pn+pd)/ms)/(r**2/gms-2*r*(mn+md))
        dmndr = 4 * np.pi * r**2 * (efpNM(pn)/ms)
        dmddr = 4 * np.pi * r**2 * (efpDM(pd)/ms)
        
        return [dpndr, dpddr,dmndr, dmddr]
    
    
    # Event to stop the solving of equations
    def event(r, y):
        if abs(y[0])>abs(y[1]):
             ev=y[0]
        else:
             ev=y[1]
        return abs(ev)-search/2
    
    def findminimum():
    
##
##

            
        
        ###################      #Some portions of the code are not included in this repository due to confidentiality.

            
        
##
##
##
        
        return frc[nmin],nmin
    
    def store(nme, dme, rnm, rdm, nres):
        # starradius = np.where(rnm >= rdm, rnm, rdm)
        # totalmass = nme+dme
        nme_frm.append(nme[nres])
        dme_frm.append(dme[nres])
        rnm_frm.append(rnm[nres])
        rdm_frm.append(rdm[nres])
        totalmass_frm.append(nme[nres]+dme[nres])
        starradius_frm.append(rdm[nres] if rdm[nres] > rnm[nres] else rnm[nres])
        deltaradius_frm.append(abs(rnm[nres]-rdm[nres]))
        
        return totalmass_frm, starradius_frm,nme_frm, rnm_frm, dme_frm, rdm_frm, deltaradius_frm
    
    # Necessary to stop the integration given the event
    event.terminal = True
    
    # range of radius of the star
    r_start = 1e-10
    r_end = 50
    step = 0.01
    search = 1e-12
    
    rlist = [[] for i in range(nstars)]
    mn = [[] for i in range(nstars)]
    md = [[] for i in range(nstars)]
    pn = [[] for i in range(nstars)]
    pd = [[] for i in range(nstars)]
    rnm = np.zeros([nstars])
    rdm = np.zeros([nstars])
    nme = np.zeros([nstars])
    dme = np.zeros([nstars])
    frc = np.zeros([nstars])
    
    def compute_tov(i, pn, pd, mn, md, nme, dme, rnm, rdm, rlist):
        init = np.array([pNM_datastars[i], pDM_datastars[i], 0, 0])
        sol = spi.solve_ivp(tov_two_fluids, [r_start, r_end], init,
                            method='RK23', max_step=step, events=event,
                            atol=1e-9, rtol=1e-9)




    
    
        ###################      #Some portions of the code are not included in this repository due to confidentiality.





    
    
        # Format the printed output
        # print(f"Iteration: {i}, RNM: {rnm[i]:.6f}")
    
    if __name__ == '__main__':
    
        num_processes = 16  # Number of parallel processes
    
        manager = multiprocessing.Manager()
        pn = manager.list([0] * nstars)
        pd = manager.list([0] * nstars)
        mn = manager.list([0] * nstars)
        md = manager.list([0] * nstars)
        rlist = manager.list([0] * nstars)
    
        nme = manager.list([0] * nstars)
        dme = manager.list([0] * nstars)
        rnm = manager.list([0] * nstars)
        rdm = manager.list([0] * nstars)
        frc = manager.list([0] * nstars)
    
    
        [frcmin,nmin]=findminimum()
        
        nsteps = 1
        cont=0
        results = []        
        ni = int(0.02*nstars)#100 #5
        nf = nmin
        nm = int(abs(nf+ni)/2)
        
        while nsteps <= niterations:
            pool = multiprocessing.Pool(processes=num_processes)
            for i in [ni, nm, nf]:
                # if frc[i] == 0:
                results.append(pool.apply_async(compute_tov, args=(
                    i, pn, pd, mn, md, nme, dme, rnm, rdm, rlist)))
    
            pool.close()
            pool.join()
     
            # Access the results
            for res in results:
                res.get()
            for i in [ni, nm, nf]:
                frc[i] = dme[i]/(nme[i]+dme[i])
    
            if frm > frc[ni] and frm > frc[nf] and nsteps == 1:
                break
            if frm < frcmin:
                break
         
            list_fr.append(round(fr, 2))
            # index_fr.append(nres) 
            
            if cont==0:     
                if frm < frc[nm]:
                    nmnew = int((nf+nm)/2)
                    ninew = nm
                    nfnew = nf
                if frm > frc[nm]:
                    nmnew = int((ni+nm)/2)
                    ninew = ni
                    nfnew = nm            
            if cont==2:     
                if frm < frc[nm]:
                    nmnew = int((ni+nm)/2)
                    ninew = ni
                    nfnew = nm
                if frm > frc[nm]:
                    nmnew = int((nm+nf)/2)
                    ninew = nm
                    nfnew = nf
                    
            print('------------')         
            print("It, FR, Radius, Mass")
            print(ni, frc[ni], rnm[ni] if rnm[ni] > rnm[ni] else rnm[ni], nme[ni]+dme[ni])
            print(nm, frc[nm], rnm[nm] if rnm[nm] > rnm[nm] else rnm[nm], nme[nm]+dme[nm])
            print(nf, frc[nf], rnm[nf] if rnm[nf] > rnm[nf] else rnm[nf], nme[nf]+dme[nf])
            print(abs(frm-frc[ni]), abs(frm-frc[nm]), abs(frm-frc[nf]))
    
            if abs(frc[nf]-frc[ni]) < search_fr:
                nres = ni
                list_fr.append(round(fr, 2))
                index_fr.append(nres)
                store(nme, dme, rnm, rdm, nres)
                cont=cont+1
                # print(nres)
    
            elif abs(frm-frc[ni]) < search_fr:
                nres = ni
                list_fr.append(round(fr, 2))
                index_fr.append(nres)
                store(nme, dme, rnm, rdm, nres)
                cont=cont+1
                # print(nres)

            elif abs(frm-frc[nm]) < search_fr:
                nres = nm
                list_fr.append(round(fr, 2))
                index_fr.append(nres)
                store(nme, dme, rnm, rdm, nres)
                cont=cont+1
                # print(nres)

            elif abs(frm-frc[nf]) < search_fr:
                nres = nf
                list_fr.append(round(fr, 2)) 
                index_fr.append(nres)
                store(nme, dme, rnm, rdm, nres)
                cont=cont+1
            else: nres=1

                # print(nres)
    
            ni = ninew
            nm = nmnew
            nf = nfnew
            print(f"Contador= {cont:.0f}")
            nsteps = nsteps+1
            
            
            if cont==1:
                cont=cont+1
                ni=nmin
                nf=nstars-1
                nm=int((nf-ni)/2)
                if frm>frc[nf]: break
                print("------------------------------------------------------------------------------")    
                
            if cont>=3: 
                break
        
        if fr<=0.96*frf:
            fr=fr+0.01
        elif fr>0.96*frf and fr<0.98*frf:
            fr=fr+0.01
        elif fr>=0.98*frf:
            fr=fr+0.01
    
        # Convert manager lists to regular lists
        pn = list(pn)
        pd = list(pd)
        mn = list(mn)
        md = list(md)
        rlist = list(rlist)
    
        nme = np.array(nme)
        dme = np.array(dme)
        rnm = np.array(rnm)
        rdm = np.array(rdm)
    
        end_time = time.time()
        elapsed_time = end_time - start_time
        print("Elapsed time:", elapsed_time)
        gc.collect()
    exportedlist = [totalmass_frm, starradius_frm,
                nme_frm, rnm_frm, dme_frm, rdm_frm, 
                deltaradius_frm]
    # Your list of data ["MT", "R", "MNM", "RNM", "MDM", "RDM", "RDM-RNM"],
    
    
    my_list = np.transpose(exportedlist)

    # Specify the file path
    frmstr = str("{:.3f}".format(frm))
    file_path = f'output.dat'
    
    #    Open the CSV file in write mode
    with open(file_path, mode='a', newline='') as file:
    # Create a CSV writer object
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(my_list)

    
    

