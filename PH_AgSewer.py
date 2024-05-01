# -*- coding: utf-8 -*-
"""V3_Calculations AGSewer.ipynb


"""

import math
import time
import numpy as np
import pandas as pd
import os
from IPython.display import HTML # Solo jupyter

#Ruta de la red de alcantarillado IMPORTANTE CAMBIAR LA RUTA POR LA CARPETA COMPARTIDA

ruta = '/content/drive/MyDrive/Colab Notebooks/AGSEWER'
#ruta = '/content/drive/MyDrive/AGSEWER'

#hola haciendo una prueba de edicion

##FUNCIONES PARA EL CALCULO HIDRAULICO

# Calculo del angulo de deflexión
def deflection_angle(x1,y1,x2,y2,x3,y3):
    return 180-math.degrees(math.acos((((x2-x1)*(x2-x3))+((y2-y1)*(y2-y3)))/(math.sqrt((x2-x1)**2+(y2-y1)**2)*math.sqrt((x2-x3)**2+(y2-y3)**2))))

##Carga de datos
#df = pd.read_csv(os.path.join(ruta,'red_alc.csv'),header=0)
def data_load():

    global df, dt, manning_cff, flow_depth, v_min, v_max, crown_exc_min, g, rh2o, frd_sbc, frd_spc, d_manhole_min, Tao_min
    #Carga de información
    df = pd.read_csv(os.path.join(ruta,'red_alc.csv'),header=0)
    dt = pd.read_csv(os.path.join(ruta,'datos_iniciales.csv'),header=0)

    #Datos iniciales
    manning_cff=float(dt['manning_cff'])
    flow_depth=float(dt['flow_depth'])
    v_min=float(dt['vel_min'])
    v_max=float(dt['vel_max'])
    crown_exc_min=float(dt['crown_exc_min'])
    g=float(dt['gravity '])
    rh2o=float(dt['density'])
    frd_sbc=float(dt['froude_sbc'])
    frd_spc=float(dt['froude_spc'])
    d_manhole_min=float(dt['d_manhole_min'])
    Tao_min=float(dt['Tao_min'])

    ##Iniciar variables
    df['x_next']=0.0
    df['y_next']=0.0
    df['Length_XY']=0.0
    df['Length_Real']=0.0
    df['abscissa']=0.0
    df['deflection_angle']=0.0
    df['S_ground']=0.0
    #Pendiente de diseño
    df['n']=manning_cff
    df['flow_depth']=flow_depth
    df['S_min']=0.0
    df['S_max']=0.0
    df['S_ref1']=0.0
    df['Theta']=0.0
    df['Area']=0.0
    df['hyd_rad']=0.0
    df['V_design']=0.0
    df['S_ref2']=0.0
    df['S_design']=0.0
    #Hydraulics
    df['Q_full']=0.0
    df['V_full']=0.0
    df['q/Q']=0.0
    df['v/V']=0.0
    df['h/D']=0.0
    df['Theta_H']=0.0
    df['Yc/D']=0.0
    df['Rh/D']=0.0
    df['V_real']=0.0
    #Altura de la lámina de agua
    df['h']=0.0
    #Profundidad Hidráulica máxima
    df['Yc']=0.0
    #Regimen de Flujo (Froude)
    df['Frd']=0.0
    df['F_R']=None
    #Altura de velocidad
    df['V^2/2g']=0.0
    #Energía específica
    df['E']=0.0
    #Fuerza Tractiva
    df['hyd_rad_des']=0.0
    df['T']=0.0
## CALCULO DE COTAS SIN PERDIDAS_______********
    df['rec_ini']=crown_exc_min
    df['rec_fin']=crown_exc_min
    #Cota Clave___________
    df['cla_ini']=0.0
    df['cla_fin']=0.0
    #Cota Batea_______
    df['bat_ini']=0.0
    df['bat_fin']=0.0
    ## PROFUNDIDAD EXACAVACIÓN COTA clave
    df['Profd_CC_ini'] = 0.0
    df['Profd_CC_fin'] = 0.0
    ## PROFUNDIDAD EXACAVACIÓN COTA BATEA
    df['Exc_ini'] = 0.0
    df['Exc_fin'] = 0.0

### COSTOSSS
    df['Media_Exc'] = 0.0
    df['Media_Exc_ft'] = 0.0


    df['Costo_tub']=0.0


    df=df.replace(np.nan, 0.0)

##Cálculos preliminares
def initializer():

    #Longitud entre pozos
    df['Length_XY']=((df['x_final']-df['x_initial'])**2+(df['y_final']-df['y_initial'])**2)**0.5
    df['abscissa']=df.groupby(['sections']).cumsum()[['Length_XY']]

    for i, row in df.iterrows():
        ## Coordenadas siguiente pozo - Requerido para el calculo del angulo de deflexion
        try:
            df.at[i,'x_next']=df.loc[np.where((df['initial_manhole']==df.at[i,'final_manhole']))]['x_final'].max()
            df.at[i,'y_next']=df.loc[np.where((df['initial_manhole']==df.at[i,'final_manhole']))]['y_final'].max()
        except:
            pass

        ## Calculo del angulo de Deflexión
        try:
            df.at[i,'deflection_angle']=deflection_angle(df.at[i,'x_initial'],df.at[i,'y_initial'],df.at[i,'x_final'],df.at[i,'y_final'],df.at[i,'x_next'],df.at[i,'y_next'])
        except:
            pass

def calculate_section(i,df_section):
    #global flow_depth, df, dt, manning_cff, flow_depth, v_min, v_max, crown_exc_min, g, rh2o, frd_sbc, frd_spc, d_manhole_min, Tao_min

    # Pendiente del terreno
    df_section.at[i,'S_ground']=(df_section.at[i,'G_Elev_Initial']-df_section.at[i,'G_Elev_Final'])/df_section.at[i,'Length_XY']
    #Pendiente mínima
    df_section.at[i,'S_min']=6.35*(v_min*manning_cff*((2*math.pi)**(2/3))*(((2*math.pi)-np.sin(2*math.pi))**(-2/3))*(df_section.at[i,'D']**(-2/3)))**2
    #Pendiente máxima
    df_section.at[i,'S_max']=6.35*(v_max*manning_cff*((2*math.pi)**(2/3))*(((2*math.pi)-np.sin(2*math.pi))**(-2/3))*(df_section.at[i,'D']**(-2/3)))**2
    #Angulo Theta para una profundidad de flujo máxima
    df_section.at[i,'Theta']=2*np.arccos(1-2*df_section.at[i,'flow_depth'])
    df_section.at[i,'Area']=((df_section.at[i,'D']**2)/8)*(df_section.at[i,'Theta']-np.sin(df_section.at[i,'Theta']))
    df_section.at[i,'hyd_rad']=(df_section.at[i,'D']/4)*(1-(np.sin(df_section.at[i,'Theta'])/df_section.at[i,'Theta']))
    df_section.at[i,'V_design']=(df_section.at[i,'Q_design']/1000)/df_section.at[i,'Area']
    df_section.at[i,'S_ref2']=6.35*(df_section.at[i,'V_design']*manning_cff*(df_section.at[i,'Theta']**(2/3))*(((df_section.at[i,'Theta'])-np.sin(df_section.at[i,'Theta']))**(-2/3))*(df_section.at[i,'D']**(-2/3)))**2


        #Pendiete segun el terreno
    if df_section.at[i,'S_min'] <= df_section.at[i,'S_ground'] and df_section.at[i,'S_max'] >= df_section.at[i,'S_ground']:
        df_section.at[i,'S_ref1'] = df_section.at[i,'S_ground']
    elif df_section.at[i,'S_min'] >= df_section.at[i,'S_ground']:
        df_section.at[i,'S_ref1'] = df_section.at[i,'S_min']
    elif df_section.at[i,'S_max'] <= df_section.at[i,'S_ground']:
        df_section.at[i,'S_ref1'] = df_section.at[i,'S_max']

        #Pendiente de Diseño

    if df_section.at[i,'S_ref2'] < df_section.at[i,'S_ref1']:
        df_section.at[i,'S_design'] = df_section.at[i,'S_ref1']
    else:
        df_section.at[i,'S_design'] = df_section.at[i,'S_ref2']


    ####Calculos Hidráulicos
    df_section.at[i,'Q_full']=(312/manning_cff)*(df_section.at[i,'D']**(8/3))*((df_section.at[i,'S_design'])**(1/2))
    df_section.at[i,'V_full']=(1/manning_cff)*((df_section.at[i,'D']/4)**(2/3))*((df_section.at[i,'S_design'])**(1/2))
    df_section.at[i,'q/Q']=(df_section.at[i,'Q_design']/1000)/(df_section.at[i,'Q_full']/1000)

        # v/V
    if df_section.at[i,'q/Q'] <= 0.06:
        df_section.at[i,'v/V'] = 10**(0.029806+0.29095*np.log10(df_section.at[i,'q/Q']))
    elif 0.06 < df_section.at[i,'q/Q'] <= 0.26:
        df_section.at[i,'v/V'] = 10**(0.013778+0.282597*np.log10(df_section.at[i,'q/Q']))
    else:
        df_section.at[i,'v/V'] = 10**(0.021763+0.289951*np.log10(df_section.at[i,'q/Q']))

        # h/D
    if df_section.at[i,'q/Q'] < 0.11:
        df_section.at[i,'h/D'] = 0.3827+0.0645*np.log(df_section.at[i,'q/Q'])
    elif 0.11 <= df_section.at[i,'q/Q'] < 0.21:
        df_section.at[i,'h/D'] = 0.60025+0.15471*np.log(df_section.at[i,'q/Q'])
    else:
        df_section.at[i,'h/D'] = 0.225+0.667*df_section.at[i,'q/Q']


    df_section.at[i,'Theta_H']=2*np.arccos(1-2*df_section.at[i,'h/D'])
    df_section.at[i,'Yc/D']=(1/8)*((df_section.at[i,'Theta_H']-np.sin(df_section.at[i,'Theta_H']))/np.sin(df_section.at[i,'Theta_H']/2))
    df_section.at[i,'Rh/D']=(1/4)*(1-(np.sin(df_section.at[i,'Theta_H'])/df_section.at[i,'Theta_H']))
    df_section.at[i,'V_real']=df_section.at[i,'v/V']*df_section.at[i,'V_full']

    #Altura de la lámina de agua
    df_section.at[i,'h']=df_section.at[i,'h/D']*df_section.at[i,'D']
    #Profundidad Hidráulica máxima
    df_section.at[i,'Yc']=df_section.at[i,'Yc/D']*df_section.at[i,'D']

    #Regimen de Flujo (Froude)
    df_section.at[i,'Frd']=df_section.at[i,'V_real']/((g*df_section.at[i,'Yc']))**0.5
    if df_section.at[i,'Frd'] <= frd_sbc:
        df_section.loc[i,'F_R']= "SUB-CRITICO"
    elif df_section.at[i,'Frd'] >= frd_spc:
        df_section.loc[i,'F_R']= "SUPER-CRITICO"
    else:
        df_section.loc[i,'F_R']= "CUASICRITICO"

    #Altura de velocidad
    df_section.at[i,'V^2/2g']=(df_section.at[i,'V_real']**2)/(2*g)
    #Energía específica
    df_section.at[i,'E']=df_section.at[i,'h']+df_section.at[i,'V^2/2g']
    #Fuerza Tractiva
    df_section.at[i,'hyd_rad_des']=df_section.at[i,'Rh/D']*df_section.at[i,'D']

    if df_section.at[i,'S_design'] <= 10:
        df_section.at[i,'T']= g*rh2o*df_section.at[i,'hyd_rad_des']*(df_section.at[i,'S_design'])
    else:
        df_section.at[i,'T']= g*rh2o*df_section.at[i,'hyd_rad_des']*np.arctan((df_section.at[i,'S_design']))

    ### CALCULO DE COTAS

    if df_section.at[i,'type'] == 1:
        df_section.at[i,'cla_ini'] = df_section.at[i,'G_Elev_Initial']-crown_exc_min
        df_section.at[i,'cla_fin'] = df_section.at[i,'cla_ini']-df_section.at[i,'Length']*(df_section.at[i,'S_design'])
    else:
        df_section.at[i,'cla_ini'] = min(df_section.loc[df_section.index[df_section['final_manhole']==df_section.at[i, 'initial_manhole']].tolist()]['cla_fin'], default=0)
        df_section.at[i,'cla_fin'] = df_section.at[i,'cla_ini']-df_section.at[i,'Length']*(df_section.at[i,'S_design'])


    df_section.at[i,'bat_ini']=df_section.at[i,'cla_ini']-df_section.at[i,'D']
    df_section.at[i,'bat_fin']=df_section.at[i,'cla_fin']-df_section.at[i,'D']

    ##PROFUNDIDAD A COTA CLAVE
    df_section.at[i,'Profd_CC_ini'] = df_section.at[i,'G_Elev_Initial'] - df_section.at[i,'cla_ini']
    df_section.at[i,'Profd_CC_fin'] = df_section.at[i,'G_Elev_Final'] - df_section.at[i,'cla_fin']
    ##PROFUNDIDAD A COTA BATEA
    df_section.at[i,'Exc_ini'] = df_section.at[i,'G_Elev_Initial'] - df_section.at[i,'bat_ini']
    df_section.at[i,'Exc_fin'] = df_section.at[i,'G_Elev_Final'] - df_section.at[i,'bat_fin']

    df_section.at[i,'Media_Exc'] = (df_section.at[i,'Exc_ini']+df_section.at[i,'Exc_fin'])/2
    df_section.at[i,'Media_Exc_ft'] = df_section.at[i,'Media_Exc']*3.28084

    ##CALCULO DE COSTOS POR TRAMO - (FORMULA COSTOS ARTICULO DEL MODELO)

    if (df_section.at[i,'D'] < 0.9144) and (df_section.at[i,'Media_Exc_ft'] < 10):
        df_section.at[i,'Costo_tub'] = (10.98*(df_section.at[i,'D']*3.28084) + 0.8*df_section.at[i,'Media_Exc_ft'] - 5.98)*(df_section.at[i,'Length_XY']*3.28084)
    elif (df_section.at[i,'D'] <= 0.9144) and (df_section.at[i,'Media_Exc_ft'] >= 10):
        df_section.at[i,'Costo_tub'] = (5.94*(df_section.at[i,'D']*3.28084) + 1.166*df_section.at[i,'Media_Exc_ft'] + 0.504*df_section.at[i,'Media_Exc_ft'] - 9.64)*(df_section.at[i,'Length_XY']*3.28084)
    elif (df_section.at[i,'D'] > 0.9144):
        df_section.at[i,'Costo_tub'] = (30.0*(df_section.at[i,'D']*3.28084) + 4.9*df_section.at[i,'Media_Exc_ft'] - 105.9)*(df_section.at[i,'Length_XY']*3.28084)

    global Total_Cost_Pipe
    Total_Cost_Pipe = np.sum(df_section['Costo_tub'])



    return df_section

#flow_depth
#list(range(70,(flow_depth)))
#flow_d = range(flow_depth)

def process():
   #itera por cada fila del dataframe
    for i, row in df.iterrows():
    

        df_s_flowdepth = df.iloc[[i]]
        df_s_flowdepth = calculate_section(i,df)

        iterations = 0

        while df_s_flowdepth.loc[i, 'V_real'] < 0.75 and df_s_flowdepth.loc[i, 'T'] < 2 or (df_s_flowdepth.loc[i, 'Frd'] > frd_sbc and df_s_flowdepth.loc[i, 'Frd'] < frd_spc):
          iterations += 1
          df.at[i,'flow_depth'] -= 0.01
          df_s_flowdepth = calculate_section(i, df)


        print (f'Tramo Tub {i+1}: Iteraciones = {iterations}')


    global manhole, Total_Cost_Manhole

    manhole = list(set(df['initial_manhole']).union(set(df['final_manhole'])))

    manhole = pd.DataFrame(manhole).drop_duplicates().rename(columns={0:"Manhole"}).reset_index(drop=True)
    manhole['max_exc'] = 0.0
    manhole['Costo_pz']=0.0
    for i, row in manhole.iterrows():
        for j, row in df.iterrows():
            if manhole.at[i,'Manhole']==df.at[j,'initial_manhole']:
                manhole.at[i,'max_exc']=max(df.at[j,'Exc_ini'],manhole.at[i,'max_exc'])
            if manhole.at[i,'Manhole']==df.at[j,'final_manhole']:
                manhole.at[i,'max_exc']=max(df.at[j,'Exc_fin'],manhole.at[i,'max_exc'])
        manhole.at[i,'Costo_pz']= 250 + (manhole.at[i,'max_exc']*3.28084)**2
    Total_Cost_Manhole = np.sum(manhole['Costo_pz'])

initial_time= time.time()
data_load()
initializer()
process()


#print("Tiempo",ruta,round((time.time()-initial_time),2), "segundos")

display(manhole)
print(f'Total Manhole Cost = {Total_Cost_Manhole}')
print(f'Total Pipe Cost = {Total_Cost_Pipe}')
print(f'Total Cost = {Total_Cost_Manhole + Total_Cost_Pipe}')


display(HTML(df.head(236).to_html(na_rep='-',formatters={'S_ground': '{:,.5%}'.format,'S_min': '{:,.5%}'.format, 'S_max': '{:,.5%}'.format,'n': '{:,.3f}'.format,'S_ref1': '{:,.5%}'.format, 'S_ref2': '{:,.5%}'.format, 'S_design': '{:,.5%}'.format})))