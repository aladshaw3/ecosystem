'''
    \file sensitivity.py
    \brief Sensitivity Analysis of the NH3 Storage and Aging Model
    \details Python script using the sensitivity.py script and object to perform a
                simple sensitivity analysis on the NH3 storage model to see which parameters
                the model is most or least sensitive to. This can be used to determine whether
                or not all the reaction schemes and aging mechanisms proposed are useful in
                determining the NH3 storage capacity on Cu-SSZ-13.
    \author Austin Ladshaw
    \date 02/13/2020
    \copyright This software was designed and built at the Oak Ridge National
                    Laboratory (ORNL) National Transportation Research Center
                    (NTRC) by Austin Ladshaw for research in the catalytic
                    reduction of NOx. Copyright (c) 2020, all rights reserved.
'''

import sensitivity as sa
import math

def NH3_Storage_Model_v0(params, conds):
    ''' The current NH3 storage model contains the following params and conds...
            params["dH1"] = reaction enthalpy (J/mol) for reaction 1
            params["dS1"] = reaction entropy (J/K/mol) for reaction 1
            params["dH2"] = reaction enthalpy (J/mol) for reaction 2
            params["dS2"] = reaction entropy (J/K/mol) for reaction 2
            params["dH3"] = reaction enthalpy (J/mol) for reaction 3
            params["dS3"] = reaction entropy (J/K/mol) for reaction 3
            params["dH4"] = reaction enthalpy (J/mol) for reaction 4
            params["dS4"] = reaction entropy (J/K/mol) for reaction 4
            params["A1"] = pre-exponential factor (1/hr/kPa^0.25) for aging reacion 1
            params["E1"] = reaction rate energy (J/mol) for aging reaction 1
            params["A2"] = pre-exponential factor (1/hr/kPa^0.25) for aging reacion 2
            params["E2"] = reaction rate energy (J/mol) for aging reaction 2
            params["A3f"] = pre-exponential factor (1/hr) for Forward aging reacion 3
            params["E3f"] = reaction rate energy (J/mol) for Forward aging reaction 3
            params["A3r"] = pre-exponential factor (1/hr) for Reverse aging reacion 3
            params["E3r"] = reaction rate energy (J/mol) for Reverse aging reaction 3
            params["Z1CuOH_o"] = initial site density (mol/L) for Z1 Cu sites
            params["Z2Cu_o"] = initial site density (mol/L) for Z2 Cu sites
            params["ZH_o"] = initial site density (mol/L) for solitary Bronsted sites
            params["ZH-ZCu_o"] = initial site density (mol/L) for Bronsted sites near inactive Z1 Cu sites
            params["ZH-CuO_o"] = initial site density (mol/L) for Bronsted sites near CuO species

            conds["T"] = temperature (K) for gas stream
            conds["P_O2"] = partial pressure (kPa) for O2 in gas stream
            conds["P_H2O"] = partial pressure (kPa) for H2O in gas stream
            conds["P_NH3"] = partial pressure (kPa) for NH3 in gas stream
            conds["aging_time"] = time spent aging (hr)
            conds["T_aging"] = temperature during aging (K) for gas stream
            conds["P_O2_aging"] = partial pressure during aging (kPa) for O2 in gas stream
            conds["P_H2O_aging"] = partial pressure during aging (kPa) for O2 in gas stream

            NOTE: aging_time = 0 for "de-greened" catalyst

            NOTE2: You can put more information in params if you want sensitivity analysis
                    to also cover the model sensitivity to things like temperature and
                    partial pressures.

            ------------ MODEL INFORMATION GIVEN BELOW -------------
            Capacity Reactions:
                (1)     (Z1CuOH) + NH3 <== ==> [(Z1CuOH)-NH3]
                (2)     (Z2Cu) + NH3 <== ==> [(Z2Cu)-NH3]
                (3)     (ZH) + NH3 <== ==> [(ZH)-NH3]
                (4)     (Z1CuOH) + H2O <== ==> [(Z1CuOH)-H2O]

            Aging Reactions:
                (1)     (ZH)(ZCu) + 0.25 O2 --> (Z2Cu) + 0.5 H2O
                (2)     (ZH) + 0.25 O2 --> (Z) + 0.5 H2O
                (3)     (Z1CuOH) <-- --> (ZH)(CuO)

            w1 = (availability of Z1CuOH after aging)
            w2 = (availability of Z2Cu after aging)
            w3 = (availability of total ZH sites after aging)
                    total ZH sites = (ZH) + (ZH)(ZCu) + (ZH)(CuO)

    '''
    R = 8.314 #J/K/mol
    #Calculate the model equilibrium parameters based on the simulation conditions given
    K1 = math.exp(-(params["dH1"]/R/conds["T"]) + (params["dS1"]/R))
    K2 = math.exp(-(params["dH2"]/R/conds["T"]) + (params["dS2"]/R))
    K3 = math.exp(-(params["dH3"]/R/conds["T"]) + (params["dS3"]/R))
    K4 = math.exp(-(params["dH4"]/R/conds["T"]) + (params["dS4"]/R))

    #Calculate the model aging parameters based on the aging conditions
    k1 = params["A1"]*math.exp(-(params["E1"]/R/conds["T_aging"]))
    k2 = params["A2"]*math.exp(-(params["E2"]/R/conds["T_aging"]))
    k3f = params["A3f"]*math.exp(-(params["E3f"]/R/conds["T_aging"]))
    k3r = params["A3r"]*math.exp(-(params["E3r"]/R/conds["T_aging"]))

    #Calculate the model site densities at the aging conditions
    A = params["Z1CuOH_o"] + params["ZH-CuO_o"]
    w1 = params["Z1CuOH_o"]*math.exp(-(k3r+k3f)*conds["aging_time"]) + ((k3r*A)/(k3r+k3f))*(1-math.exp(-(k3r+k3f)*conds["aging_time"]))
    w2 = params["Z2Cu_o"] + params["ZH-ZCu_o"]*(1-math.exp(-k1*math.pow(conds["P_O2_aging"],0.25)*conds["aging_time"]))
    w3 = params["ZH-ZCu_o"]*math.exp(-k1*math.pow(conds["P_O2_aging"],0.25)*conds["aging_time"])
    w3 += params["ZH_o"]*math.exp(-k2*math.pow(conds["P_O2_aging"],0.25)*conds["aging_time"])
    w3 += (A - w1)

    #Calculate the storage capacity
    NH3 = w1*(K1*conds["P_NH3"])/(1+K1*conds["P_NH3"]+K4*conds["P_H2O"])
    NH3 += w2*(K2*conds["P_NH3"])/(1+K2*conds["P_NH3"])
    NH3 += w3*(K3*conds["P_NH3"])/(1+K3*conds["P_NH3"])

    return NH3
#End NH3_Storage_Model_v0

def NH3_Storage_Model_v0_1(params, conds):
    ''' The next NH3 storage model contains the following params and conds...
            params["dH1"] = reaction enthalpy (J/mol) for reaction 1
            params["dS1"] = reaction entropy (J/K/mol) for reaction 1
            params["dH2"] = reaction enthalpy (J/mol) for reaction 2
            params["dS2"] = reaction entropy (J/K/mol) for reaction 2
            params["dH3"] = reaction enthalpy (J/mol) for reaction 3
            params["dS3"] = reaction entropy (J/K/mol) for reaction 3
            params["dH4"] = reaction enthalpy (J/mol) for reaction 4
            params["dS4"] = reaction entropy (J/K/mol) for reaction 4
            params["k1"] = aging rate (1/hr/kPa^0.25) for aging reacion 1
            params["k2"] = aging rate (1/hr/kPa^0.25) for aging reacion 2
            params["k3f"] = aging rate (1/hr) for Forward aging reacion 3
            params["k3r"] = aging rate (1/hr) for Reverse aging reacion 3
            params["Z1CuOH_o"] = initial site density (mol/L) for Z1 Cu sites
            params["Z2Cu_o"] = initial site density (mol/L) for Z2 Cu sites
            params["ZH_o"] = initial site density (mol/L) for solitary Bronsted sites
            params["ZH-ZCu_o"] = initial site density (mol/L) for Bronsted sites near inactive Z1 Cu sites
            params["ZH-CuO_o"] = initial site density (mol/L) for Bronsted sites near CuO species

            conds["T"] = temperature (K) for gas stream
            conds["P_O2"] = partial pressure (kPa) for O2 in gas stream
            conds["P_H2O"] = partial pressure (kPa) for H2O in gas stream
            conds["P_NH3"] = partial pressure (kPa) for NH3 in gas stream
            conds["aging_time"] = time spent aging (hr)
            conds["T_aging"] = temperature during aging (K) for gas stream
            conds["P_O2_aging"] = partial pressure during aging (kPa) for O2 in gas stream
            conds["P_H2O_aging"] = partial pressure during aging (kPa) for O2 in gas stream

            NOTE: aging_time = 0 for "de-greened" catalyst

            NOTE2: You can put more information in params if you want sensitivity analysis
                    to also cover the model sensitivity to things like temperature and
                    partial pressures.

            ------------ MODEL INFORMATION GIVEN BELOW -------------
            Capacity Reactions:
                (1)     (Z1CuOH) + NH3 <== ==> [(Z1CuOH)-NH3]
                (2)     (Z2Cu) + NH3 <== ==> [(Z2Cu)-NH3]
                (3)     (ZH) + NH3 <== ==> [(ZH)-NH3]
                (4)     (Z1CuOH) + H2O <== ==> [(Z1CuOH)-H2O]

            Aging Reactions:
                (1)     (ZH)(ZCu) + 0.25 O2 --> (Z2Cu) + 0.5 H2O
                (2)     (ZH) + 0.25 O2 --> (Z) + 0.5 H2O
                (3)     (Z1CuOH) <-- --> (ZH)(CuO)

            w1 = (availability of Z1CuOH after aging)
            w2 = (availability of Z2Cu after aging)
            w3 = (availability of total ZH sites after aging)
                    total ZH sites = (ZH) + (ZH)(ZCu) + (ZH)(CuO)

    '''
    R = 8.314 #J/K/mol
    #Calculate the model equilibrium parameters based on the simulation conditions given
    K1 = math.exp(-(params["dH1"]/R/conds["T"]) + (params["dS1"]/R))
    K2 = math.exp(-(params["dH2"]/R/conds["T"]) + (params["dS2"]/R))
    K3 = math.exp(-(params["dH3"]/R/conds["T"]) + (params["dS3"]/R))
    K4 = math.exp(-(params["dH4"]/R/conds["T"]) + (params["dS4"]/R))

    #Calculate the model aging parameters based on the aging conditions
    k1 = params["k1"]
    k2 = params["k2"]
    k3f = params["k3f"]
    k3r = params["k3r"]

    #Calculate the model site densities at the aging conditions
    A = params["Z1CuOH_o"] + params["ZH-CuO_o"]
    w1 = params["Z1CuOH_o"]*math.exp(-(k3r+k3f)*conds["aging_time"]) + ((k3r*A)/(k3r+k3f))*(1-math.exp(-(k3r+k3f)*conds["aging_time"]))
    w2 = params["Z2Cu_o"] + params["ZH-ZCu_o"]*(1-math.exp(-k1*math.pow(conds["P_O2_aging"],0.25)*conds["aging_time"]))
    w3 = params["ZH-ZCu_o"]*math.exp(-k1*math.pow(conds["P_O2_aging"],0.25)*conds["aging_time"])
    w3 += params["ZH_o"]*math.exp(-k2*math.pow(conds["P_O2_aging"],0.25)*conds["aging_time"])
    w3 += (A - w1)

    #Calculate the storage capacity
    NH3 = w1*(K1*conds["P_NH3"])/(1+K1*conds["P_NH3"]+K4*conds["P_H2O"])
    NH3 += w2*(K2*conds["P_NH3"])/(1+K2*conds["P_NH3"])
    NH3 += w3*(K3*conds["P_NH3"])/(1+K3*conds["P_NH3"])

    return NH3
#End NH3_Storage_Model_v0_1

# ------------- Run Tests and Simulations --------------
params = {}
conds_lb = {}
conds_ub = {}

params["dH1"] = -87325.61
params["dS1"] = -6.1568
params["dH2"] = -79847.49
params["dS2"] = -117.1946
params["dH3"] = -88938.66
params["dS3"] = -101.1138
params["dH4"] = -34673.53
params["dS4"] = 79.5312

#Set for NH3_Storage_Model_v0
#params["A1"] = 1727.17
#params["E1"] = 70409.04
#params["A2"] = 26.457
#params["E2"] = 48330.623
#params["A3f"] = 0.097438
#params["E3f"] = 1595.34
#params["A3r"] = 0.1640
#params["E3r"] = -4530.22

#Set for NH3_Storage_Model_v0_1
params["k1"] = 0.6458224
params["k2"] = 0.1174927
params["k3f"] = 0.0814843
params["k3r"] = 0.2725287

params["Z1CuOH_o"] = 0.050158
params["Z2Cu_o"] = 0.0355518
params["ZH_o"] = 0.0140659
params["ZH-ZCu_o"] = 0.0166239
params["ZH-CuO_o"] = 0.0087223

conds_lb["T"] = 423.15
conds_lb["P_O2"] = 0.188292810335608
conds_lb["P_H2O"] = 2.73842147591615
conds_lb["P_NH3"] = 0.001271108296667
conds_lb["aging_time"] = 0.0
conds_lb["T_aging"] = 1073.15
conds_lb["P_O2_aging"] = 9.415
conds_lb["P_H2O_aging"] = 9.415

conds_ub["T"] = 623.15
conds_ub["P_O2"] = 0.188292810335608
conds_ub["P_H2O"] = 10.5618600586355
conds_ub["P_NH3"] = 0.09811135049
conds_ub["aging_time"] = 16.0
conds_ub["T_aging"] = 1073.15
conds_ub["P_O2_aging"] = 9.415
conds_ub["P_H2O_aging"] = 9.415

conds_tuples = {}
for item in conds_lb:
    conds_tuples[item] = (conds_lb[item], conds_ub[item])

#analysis = sa.SensitivitySweep(NH3_Storage_Model_v0, params, conds_tuples)
analysis = sa.SensitivitySweep(NH3_Storage_Model_v0_1, params, conds_tuples)
file_name_simple = "NH3-Analysis-Results-Simple.txt"
file_name_full = "NH3-Analysis-Results-Exhaustive.txt"
rel = True
per = 10
#analysis.run_sweep(file_name_simple, rel, per)
analysis.run_exhaustive_sweep(file_name_full, rel, per)

#print(analysis)
