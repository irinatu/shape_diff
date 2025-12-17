#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 00:03:53 2023

@author: ituszynska

Compare two shape files in .map formats and produce pictures with smooth reactivities and thair correlation

"""
import argparse, sys, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


#Definition of parameters and arguments
def pars_inp():
    parser = argparse.ArgumentParser(
        description=(
            "Compare two SHAPE reactivity .map files (in vitro vs in cellulo), "
            "smooth reactivities, calculate differences, correlations, and "
            "generate plots and VARNA visualizations."
        )
    )

    parser.add_argument(
        "-d", "--dot",
        dest="Dot",
        type=str,
        required=True,
        help="Secondary structure file in dot-bracket format"
    )

    parser.add_argument(
        "-v", "--vitro",
        dest="Vitro",
        type=str,
        required=True,
        help="In vitro SHAPE reactivity file (.map)"
    )

    parser.add_argument(
        "-c", "--cell",
        dest="Cell",
        type=str,
        required=True,
        help="In cellulo SHAPE reactivity file (.map)"
    )

    parser.add_argument(
        "-t", "--vitro-name",
        dest="v_name",
        type=str,
        default="",
        help="Label for in vitro SHAPE data (used in plots/output filenames)"
    )

    parser.add_argument(
        "-e", "--cell-name",
        dest="c_name",
        type=str,
        default="",
        help="Label for in cellulo SHAPE data (used in plots/output filenames)"
    )

    parser.add_argument(
        "-o", "--out",
        dest="Out",
        action="store_true",
        help="Save tab-delimited output with calculated statistics"
    )

    parser.add_argument(
        "-s", "--shape",
        dest="Shape",
        action="store_true",
        help="Save significant SHAPE differences in VARNA-compatible format"
    )

    parser.add_argument(
        "-f", "--figures",
        dest="Figures",
        action="store_true",
        help=(
            "Generate VARNA structure figures (requires varnaapi and VARNA .jar). "
            "Recommended for RNAs up to ~200 nt."
        )
    )

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()

#loocking for the VARNA path from enviroment variable VARNA_JAR
def get_varna_path():
    varna = os.environ.get("VARNA_JAR")
    if varna is None:
        sys.exit(
            "ERROR: VARNA_JAR environment variable not set.\n"
            "Please set it, e.g.:\n"
            "  export VARNA_JAR=/path/to/VARNAv3-94a.jar"
        )
    if not os.path.isfile(varna):
        sys.exit(f"ERROR: VARNA_JAR points to non-existing file: {varna}")
    return varna

#checking of dot-bracket file
def check(lista, expected):
    unexpected = []
    
    for char in lista:
        if char not in expected:
            unexpected.append(char)
    
    if unexpected:
        print("warning: Unexpected signs in your dot bracket file, expected are ", expected)
        for znak in unexpected:
            print(znak)
    else:
        pass

# significans of reactivities changes 
def significantly_changed(df):
    if (df['Z'] > 0.0) & (df['S'] >= 1.0): return df['dReac'].values
    else: return 0.0 

#fill with zeros nan values
def fillZero(tab):
    tab['smoth_reac_plot'] = tab['smoth_reac']
    tab = tab.apply(lambda col: col.fillna(0.0) if col.name not in ['dReac', 'smoth_reac_plot'] else col)
    return tab

opts= pars_inp()

#setting up the name of output files
if opts.v_name == '': vitro = os.path.splitext(os.path.basename(opts.Vitro))[0]
else: vitro = opts.v_name

if opts.c_name == '': cell = os.path.splitext(os.path.basename(opts.Cell))[0]
else: cell = opts.c_name
    
#open and smooth reactivities from the first .map file
vit = pd.read_csv(opts.Vitro, delimiter="\t", header = None)
vit.columns = ["nr", "react", "SE", "nt"]
vit['react'] = np.where(vit['react'] == -999, np.nan, vit['react'])
vit['SE'] = np.where(vit['react'].isna(), np.nan, vit['SE'])
# Converting 'MaskedArray' to Pandas Series
react_series = pd.Series(np.ma.masked_invalid(vit['react']).data)

# Applying rolling mean to the reactivities
vit['smoth_reac'] = react_series.rolling(3, min_periods=1, center=True).mean()
vit['smoth_reac'] = np.where(vit['react'].isna(), np.nan, vit['smoth_reac'])

react_series_s = pd.Series(np.ma.masked_invalid(vit['SE']).data)
vit['smoth_SE'] = (np.sqrt(react_series_s.rolling(3, min_periods=1, center=True).apply(lambda x: np.sum(x**2))))/3


#open and smooth reactivities from the second .map file
cel = pd.read_csv(opts.Cell, delimiter="\t", header = None)
cel.columns = ["nr", "react", "SE", "nt"]
#convert -999 to NAN
cel['react'] = np.where(cel['react'] == -999, np.nan, cel['react'])
cel['SE'] = np.where(cel['react'].isna(), np.nan, cel['SE'])

#Mask invalid series
react_series_c = pd.Series(np.ma.masked_invalid(cel['react']).data)

# Applying rolling mean to the Series
cel['smoth_reac'] = react_series_c.rolling(3, min_periods=1, center=True).mean()
cel['smoth_reac'] = np.where(cel['react'].isna(), np.nan, cel['smoth_reac'])

#mask and calculate SE_smooth like 1/3sqrt(spred**2+s**2+spo**2)
react_series_sc = pd.Series(np.ma.masked_invalid(cel['SE']).data)

cel['smoth_SE'] = (np.sqrt(react_series_sc.rolling(3, min_periods=1, center=True).apply(lambda x: np.sum(x**2))))/3


#calculation of the statistics

cel['dReac'] = vit['smoth_reac'] - cel['smoth_reac'] 
#print ("REACT", cel["dReac"][200], vit['smoth_reac'][200], cel['smoth_reac'][200], vit['react'][200],  cel['react'][200])
cel['Z'] = 1-(1.96*(vit['smoth_SE']+cel['smoth_SE'])/abs(cel['dReac']))


mean_dR = cel['dReac'].mean()
std_dR = cel['dReac'].std()
cel['S'] = (cel['dReac'] - mean_dR)/std_dR

#co najmnien 3 z okienka 5 maja |S|>=1 i Z > 0.0
cel.loc[((cel['Z'] > 0.0) & (abs(cel['S']) >= 1.0)).rolling(5, center = True).sum() >=3,'dif_sig'] = cel['dReac']
cel['dif_sig'].fillna(0.0, inplace=True)


# Tworzenie nowej tablicy
corel_table = pd.DataFrame({'cel': cel['smoth_reac'], 'vit': vit['smoth_reac']})

# Usuwanie wierszy zawierających przynajmniej jedną wartość NaN
corel_table = corel_table.dropna()

# Obliczanie korelacji pomiędzy kolumnami 'cel' i 'vit'
r_s, p_s = spearmanr(corel_table['cel'], corel_table['vit'])


#save txt output with calculated statistics
if opts.Out:
    data = {
    'Nr': vit['nr'],
    'Seq': vit ['nt'],
    'deltaShape': cel['dif_sig'].round(2),
    'Z-factor': cel['Z'].round(2),	
    'Std_Score': cel['S'].round(2),
    vitro + '_react': vit['react'],
    vitro + '_smooth_rea': vit['smoth_reac'],
    cell + '_react': cel['react'],
    cell + '_smooth_rea': cel['smoth_reac']
    }

    df_n = pd.DataFrame(data)
    df_n = df_n.fillna(-999)

    # Zapis do pliku CSV
    df_n.to_csv(vitro + '_' + cell + '.out', index=False, sep='\t')


vit = fillZero(vit)
cel = fillZero(cel)

if opts.Shape:
    data = {
    'Nr': vit['nr'],
    'deltaShape': cel['dif_sig'].round(2),
    }

    df_s = pd.DataFrame(data)

    # Zapis do pliku CSV
    df_s.to_csv(vitro + '_' + cell + 'diff.shape', index=False, sep='\t', header=False)

if opts.Figures:
    import varnaapi

    #pictures of the secondary structure with shape reactivities 
    varnaapi.set_VARNA(get_varna_path())
    
    
    color_shape_v = {vit['react'].min():'#000000',0.40:'#000000',0.41:'#FFCC00',0.80:'#FFCC00',0.81:'#990000', vit['react'].max():'#990000'}
    color_shape_c = {cel['react'].min():'#000000',0.40:'#000000',0.41:'#FFCC00',0.80:'#FFCC00',0.81:'#990000', cel['react'].max():'#990000'}
    color_shape_v_sm = {vit['smoth_reac'].min():'#000000',0.40:'#000000',0.41:'#FFCC00',0.80:'#FFCC00',0.81:'#990000', vit['smoth_reac'].max():'#990000'}
    color_shape_c_sm = {cel['smoth_reac'].min():'#000000',0.40:'#000000',0.41:'#FFCC00',0.80:'#FFCC00',0.81:'#990000', cel['smoth_reac'].max():'#990000'}
    #color_shape_d = {cel['dif_sig'].min():'#990033', -0.03:'#FFFFFF', 0.03:'#FFFFFF', cel['dif_sig'].max():'#000099'}
    color_shape_d = {-1:'#990033', 0:'#FFFFFF', 1:'#000099'}
    
    #open file with secondary structure in dot bracket format
    with open(opts.Dot, "r") as file:
        lines = file.readlines()
    
    # Odczytanie drugiej linii
    sequence = lines[1].strip()
    check(sequence, ['A', 'C', 'G', 'T', 'U'])
    
    # Odczytanie trzeciej linii
    structure = lines[2].strip()
    check(structure, ['.', '(', ')'])
    
    
    for sh, sh_c, n in zip((vit['react'], cel['react'], vit['smoth_reac'], cel['smoth_reac'], cel['dif_sig']), 
                           (color_shape_v, color_shape_c, color_shape_v_sm, color_shape_c_sm, color_shape_d),
                           (vitro, cell, vitro + "_smoth", cell + "_smoth", vitro + '_' + cell + "_"+"diff")):
        #print ("UWAGA",n, sh, cel['dif_sig'].max(), max(sh))
        v = varnaapi.Structure(sequence ,  structure)
        if n ==  vitro + '_' + cell + "_"+"diff":
            tres = max(abs(min(sh)), abs(max(sh)))
            print ("UWAGA", tres,min(sh), max(sh) )
            v.add_colormap(sh, vMin = -tres, vMax = tres, style = sh_c)
        else: v.add_colormap(sh, vMin = min(sh), vMax = max(sh), style = sh_c)
        v.update( resolution = 4, spaceBetweenBases = 0.7)
        v.savefig( n +'_shape' +'.png')
    

fig, axes = plt.subplots(nrows=3, ncols=1)
fig.tight_layout()
font_size = 5

cel['X'] = list(range(1, len(cel['S'])+1))
vit['X'] = list(range(1, len(vit['smoth_reac'])+1))

ax1 = vit.plot(x = 'X', y = 'smoth_reac_plot', drawstyle='steps-mid', color='blue', ax=axes[0], xlabel = '', ylabel = 'reactivities', fontsize = font_size)
ax1 =cel.plot(x = 'X', y = 'smoth_reac_plot', drawstyle='steps-mid', color='red', ax=axes[0], xlabel = '', fontsize = font_size)
leg =ax1.legend([vitro, cell], fontsize=font_size)
ax1.tick_params(axis='both', labelsize=font_size)
ax1.yaxis.label.set_size(font_size)

ax2 = cel.plot(x = 'X', y = 'dReac', drawstyle='steps-mid', color='black', ax=axes[1], xlabel = 'nr of nucleotides',ylabel = 'react. diff', fontsize=font_size)
leg2 = ax2.legend([vitro + " - " + cell], fontsize=font_size)
#leg2.prop.set_size(font_size)
ax2.fill_between(cel['X'], cel['dif_sig'], where=(cel['dif_sig'] >= 0.0), color='green')
ax2.fill_between(cel['X'], cel['dif_sig'], where=(cel['dif_sig'] < 0.0), color='blue')
ax2.axhline(0, color='black')
ax2.tick_params(axis='both', labelsize=font_size)
ax2.xaxis.label.set_size(font_size)
ax2.yaxis.label.set_size(font_size)

# the same scale on X for plot nr 1 and 2
axes[1].set_xlim(axes[0].get_xlim())

ax3 = axes[2]
ax3.plot(vit['smoth_reac'], cel['smoth_reac'], 's', color='red', markersize=2) 
ax3.set_xlabel(vitro, fontsize=font_size)
ax3.set_ylabel(cell, fontsize=font_size)
leg3 = ax3.legend ([ f'R_Sp={r_s:.2f}, p_Sp={p_s:.2f}' ], fontsize=font_size)
#leg3.prop.set_size(font_size)
ax3.tick_params(axis='both', labelsize=font_size)
plt.rcParams.update({'font.size': 3})
plt.tight_layout()

#plt.show()
plt.savefig(vitro + "_"+ cell + "_comp.png", dpi=1200)

