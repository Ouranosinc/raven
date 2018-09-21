import numpy as np

def simulation(data, params):
    '''
    Cema-Neige snow model

    Input:
    1. Data - pandas dataframe with correspondent timeseries:
        'Temp'- mean daily temperature (Celsium degrees)
        'Prec'- mean daily precipitation (mm/day)
    2. Params - list of model parameters:
        'CTG' - dimensionless weighting coefficient of the snow pack thermal state
                [0, 1]
        'Kf'  - day-degree rate of melting (mm/(day*celsium degree))
                [1, 10]
    Output:
    Total amount of liquid and melting precipitation daily timeseries
    (for coupling with hydrological model)
    '''

    ### reading the data ###
    Temp = data['Temp']
    Prec = data['Prec']
    FraqSolidPrecip = np.where(Temp < -0.2, 1, 0)

    CTG, Kf = params

    ### initialization ###
    ## constants ##
    # melting temperature
    Tmelt = 0
    # Threshold for solid precip
    # function for Mean Annual Solid Precipitation
    def MeanAnnualSolidPrecip(data):
        annual_vals = [data.Prec.ix[data.Prec.index.year == i][data.Temp < -0.2].sum()\
                       for i in np.unique(data.index.year)]
        return np.mean(annual_vals)

    MASP = MeanAnnualSolidPrecip(data)
    Gthreshold = 0.9*MASP
    MinSpeed = 0.1

    ## model states ##
    G = 0
    eTG = 0
    PliqAndMelt = 0

    ### ouput of snow model
    PliqAndMelt = np.zeros(len(Temp))

    for t in range(len(Temp)):
        ### solid and liquid precipitation accounting
        # liquid precipitation
        Pliq = (1 - FraqSolidPrecip[t]) * Prec[t]
        # solid precipitation
        Psol = FraqSolidPrecip[t] * Prec[t]
        ### Snow pack volume before melt
        G = G + Psol
        ### Snow pack thermal state before melt
        eTG = CTG * eTG + (1 - CTG) * Temp[t]
        # control eTG
        if eTG > 0: eTG = 0
        ### potential melt
        if (int(eTG) == 0) & (Temp[t] > Tmelt):
            PotMelt = Kf * (Temp[t] - Tmelt)
            if PotMelt > G: PotMelt = G
        else:
            PotMelt = 0
        ### ratio of snow pack cover (Gratio)
        if G < Gthreshold:
            Gratio = G/Gthreshold
        else:
            Gratio = 1
        ### actual melt
        Melt = ((1 - MinSpeed) * Gratio + MinSpeed) * PotMelt
        ### snow pack volume update
        G = G - Melt
        ### Gratio update
        if G < Gthreshold:
            Gratio = G/Gthreshold
        else:
            Gratio = 1

        ### Water volume to pass to the hydrological model
        PliqAndMelt[t] = Pliq + Melt

    return PliqAndMelt

def bounds():
    '''
    'CTG' - dimensionless weighting coefficient of the snow pack thermal state
            [0, 1]
    'Kf'  - day-degree rate of melting (mm/(day*celsium degree))
            [1, 10]
    '''
    bnds = ((0, 1), (1, 10))
    return bnds
