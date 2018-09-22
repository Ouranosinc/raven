import numpy as np
import cema_neige

def simulation(data, params):
    '''
    SIMHYD hydrological model (Chiew, 2009)
    coupled with Cema-Neige snow model

    Input:
    data - pandas dataframe with correspondent variables:
        'Temp' - daily temperature (Celsium degrees)
        'Prec' - daily precipitation (mm/day)
        'Evap' - potential evapotraspiration (mm/day)
    parameters:
        'INSC' - interception store capacity (mm)
            [0, 50]
        'COEFF'- maximum infiltration loss
            [0, 400]
        'SQ'   - infiltration loss exponent
            [0, 10]
        'SMSC' - soil moisture storage capacity
            [1, 1000]
        'SUB'  - constant of proportionality in interflow equation
            [0, 1]
        'CRAK' - constant of proportionality in groundwater recharge equation
            [0, 1]
        'K'    - baseflow linear recession parameter
            [0, 1]
        'etmul'- added parameter to convert maxT to PET
            [0.1, 3]
        ### Muskinghum routing parameters
        'DELAY'- runoff delay
            [0.1, 5]
        'X_m'  - transformation parameter
            [0.01, 0.5]
        ### Cema-Neige parameters
        X5 : dimensionless weighting coefficient of the snow pack thermal state
            [0, 1]
        X6 : day-degree rate of melting (mm/(day*celsium degree))
            [1, 10]
    '''
    ### read parameters ###
    INSC, COEFF, SQ, SMSC, SUB, CRAK, K, etmul, DELAY, X_m, X5, X6 = params

    ### read the data ###
    #Temp = data['Temp']
    #Prec = data['Prec']
    Evap = data['Evap'] * etmul
    Prec = cema_neige.simulation(data, [X5, X6])

    ### states and parameters initialization ###
    # total runoff
    U    = np.zeros(len(Prec))
    # interception store
    IMAX = np.zeros(len(Prec))
    # interception amount
    INT  = np.zeros(len(Prec))
    # interception runoff
    INR  = np.zeros(len(Prec))
    # infiltration capacity
    RMO  = np.zeros(len(Prec))
    # direct runoff
    IRUN = np.zeros(len(Prec))
    # Soil evaporation
    ET   = np.zeros(len(Prec))
    # Saturation excess runoff and interflow
    SRUN = np.zeros(len(Prec))
    # Recharge
    REC  = np.zeros(len(Prec))
    # Infiltration into soil store
    SMF  = np.zeros(len(Prec))
    # potential evapotranspiration (PET - interception)
    POT  = np.zeros(len(Prec))
    # baseflow
    BAS  = np.zeros(len(Prec))
    # soil moisture storage
    SMS  = np.zeros(len(Prec))
    # ground water storage
    GW   = np.zeros(len(Prec))

    GWt1, GWt0 = 0, 0
    SMSt0 = 0.5
    SMSt1 = SMSt0 * SMSC

    for t in range(len(Prec)):
        # calculate interception store
        IMAX[t] = min(INSC, Evap[t])
        # then calculate interception
        INT[t] = min(IMAX[t], Prec[t])
        # calculate runoff after interception
        INR[t] = Prec[t] - INT[t]
        # calculate infiltration capacity
        RMO[t] = min(COEFF*np.exp(-SQ*SMSt1/SMSC), INR[t])
        # calculate direct runoff after loading to infiltration capacity
        IRUN[t] = INR[t] - RMO[t]
        # saturation excess runoff and interflow
        SRUN[t] = SUB * SMSt1 / SMSC * RMO[t]
        # calculate recharge
        REC[t] = CRAK * SMSt1 / SMSC * (RMO[t] - SRUN[t])
        # calculate infiltration into soil store
        SMF[t] = RMO[t] - SRUN[t] - REC[t]
        # calculate potential ET (amount of Evap after loses)
        POT[t] = Evap[t] - INT[t]
        # calculate soil evaporation
        ET[t] = min(10 * SMSt1/SMSC, POT[t])
        # calculate soil moisture storage (SMS) overflow
        SMS[t] = SMSt1 + SMF[t] - ET[t]
        # update states of SMS, REC and SMSt1
        if SMS[t] > SMSC:
            SMS[t] = SMSC
            REC[t] = REC[t] + SMS[t] - SMSC
        SMSt1 = SMS[t]
        # calculate baseflow
        BAS[t] = K * GWt1
        # calculate ground water storage
        GW[t] = GWt1 + REC[t] - BAS[t]
        # update state of GWt1
        GWt1 = GW[t]
        # final runoff (effective precipitation) calculation
        U[t] = IRUN[t] + SRUN[t] + BAS[t]

    ### Muskinghum routing scheme ###
    # initialize transformed runoff
    Q = np.zeros(len(U))
    # calculate Muskinghum components
    if (2*DELAY*X_m < 1) & (2*DELAY*(1-X_m) > 1):
        C0 = (-DELAY*X_m+0.5)/(DELAY*(1-X_m)+0.5)
        C1 = (DELAY*X_m+0.5)/(DELAY*(1-X_m)+0.5)
        C2 = (DELAY*(1-X_m)-0.5)/(DELAY*(1-X_m)+0.5)
    else:
        C0 = 0
        C1 = 1
        C2 = 0
    # check formal relations
    if (C0 + C1 + C2) != 1.0:
        C0 = 0
        C1 = 1
        C2 = 0
    # start transformation
    Q[0] = U[0]
    for t in range(len(U)-1):
        Q[t+1] = C0 * U[t+1] + C1 * U[t] + C2 * Q[t]
        # control Q
        if Q[t+1] < 0: Q[t+1] = 0

    return Q

def bounds():
    '''
    'INSC' - interception store capacity (mm)
        [0, 50]
    'COEFF'- maximum infiltration loss
        [0, 400]
    'SQ'   - infiltration loss exponent
        [0, 10]
    'SMSC' - soil moisture storage capacity
        [1, 1000]
    'SUB'  - constant of proportionality in interflow equation
        [0, 1]
    'CRAK' - constant of proportionality in groundwater recharge equation
        [0, 1]
    'K'    - baseflow linear recession parameter
        [0, 1]
    'etmul'- added parameter to convert maxT to PET
        [0.1, 3]
    ### Muskinghum routing parameters
    'DELAY'- runoff delay
        [0.1, 5]
    'X_m'  - transformation parameter
        [0.01, 0.5]
    ### Cema-Neige parameters
    X5 : dimensionless weighting coefficient of the snow pack thermal state
        [0, 1]
    X6 : day-degree rate of melting (mm/(day*celsium degree))
        [1, 10]
    '''
    bnds = ((0, 50), (0, 400), (0, 10), (0, 1000),\
            (0, 1), (0, 1), (0, 1), (0.1, 3),\
            (0.1, 5), (0.01, 0.5), (0, 1), (1, 10))
    return bnds

# import modules for interaction()
import pandas as pd
import sys
sys.path.append('../tools/')
from wfdei_to_lumped_dataframe import dataframe_construction
from metrics import NS

def interaction(river_name, path_to_scheme, path_to_observations,\
    INSC, COEFF, SQ, SMSC, SUB, CRAK, K, etmul, DELAY, X_m, X5, X6):

    # simulate our modeled hydrograph
    data = dataframe_construction(path_to_scheme)
    data['Qsim'] = simulation(data, [INSC, COEFF, SQ, SMSC, SUB, CRAK, K,\
    etmul, DELAY, X_m, X5, X6])

    # read observations
    obs = pd.read_csv(path_to_observations, index_col=0, parse_dates=True,
                      squeeze=True, header=None, names=['Date', 'Qobs'])

    # concatenate data
    data = pd.concat([data, obs], axis=1)

    # calculate efficiency criterion
    # slice data only for observational period and drop NA values
    data_for_obs = data.ix[obs.index, ['Qsim', 'Qobs']].dropna()
    eff = NS(data_for_obs['Qobs'], data_for_obs['Qsim'])

    # plot
    ax = data.ix[obs.index, ['Qsim', 'Qobs']].plot(figsize=(10, 7), style=['b-', 'k.'])
    ax.set_title(river_name + ' daily runoff modelling, ' + 'Nash-Sutcliffe efficiency: {}'.format(np.round(eff, 2)))
