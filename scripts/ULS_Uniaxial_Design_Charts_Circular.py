# ==============================================================================
# Copyright (c) 2026 I. Mavropodis, Ch. Zeris. All rights reserved.
# 
# This source code is part of the Diploma Thesis developed at the 
# National Technical University of Athens (NTUA).
# It is provided for review and academic purposes only. 
# Unauthorized copying, modification, or commercial use is strictly prohibited.
# ==============================================================================
import numpy
data = {}
eps2 = 2.0         #1/1000
ecu2 = 3.5         #1/1000
Es = 200000        #MPa
fyk = 500          #MPa
eyk = fyk/Es*1000  #1/1000
gammas = 1.15
eyd = eyk/gammas   #1/1000
eud = 200          #1/1000
def concrete(strains):
    sc_fcd = numpy.zeros_like(strains)
    parabola = (strains>0) & (strains<=eps2)
    sc_fcd[parabola] = 1-(1-strains[parabola]/eps2)**2
    rectangular = (strains>eps2)
    sc_fcd[rectangular] = 1
    return sc_fcd
def rebar(strains):
    ss_fyd = numpy.zeros_like(strains)
    elastic = numpy.abs(strains)<=eyd
    ss_fyd[elastic] = strains[elastic]/eyd
    plastic = (numpy.abs(strains)>eyd)
    ss_fyd[plastic] = numpy.sign(strains[plastic])
    return ss_fyd
d1_h_list = [x/100 for x in range(5,21,5)]
omegatot_list = numpy.arange(0,2.05,0.05)
bars = 8
for d1_h in d1_h_list:
    data[d1_h] = {'curves': {}, 'lines': {}}
    ec1_list = []
    ec2_list = []
    step = 0.001
    multiplier = int(1/step)
    r_start_1 = 0
    r_end_1 = (-eud-ecu2)/(1-d1_h)
    start_1 = int(round(r_start_1*multiplier))
    finish_1 = int(round(r_end_1*multiplier))
    for r in range(start_1,finish_1-1,-1):
        r = r/float(multiplier)
        es1i = -eud
        ec1i = es1i+r*d1_h
        ec2i = es1i-r*(1-d1_h)
        ec1_list.append(ec1i)
        ec2_list.append(ec2i)
    r_start_2 = (-eud-ecu2)/(1-d1_h)
    r_end_2 = 0
    start_2 = int(round(r_start_2*multiplier))
    finish_2 = int(round(r_end_2*multiplier))
    for r in range(start_2,finish_2+1,1):
        r = r/float(multiplier)
        ec2i = ecu2
        ec1i = ec2i+r
        ec1_list.append(ec1i)
        ec2_list.append(ec2i)
    grid = 200
    rsi = 0.5-d1_h
    angles = numpy.linspace(0,2*numpy.pi,bars,endpoint=False)
    ysi = 0.5+rsi*numpy.sin(angles)
    dy = 1/grid
    yci = numpy.linspace(dy/2,1-dy/2,grid)
    bci = 2*numpy.sqrt(0.5**2-(yci-0.5)**2)
    Ac = numpy.pi*(0.5**2)
    dAc = (bci*dy)/Ac
    for omegatot in omegatot_list:    
        vd_list = []
        md_list = []
        omegai = numpy.full(bars,omegatot/bars)
        ec_amount = range(len(ec2_list))
        for i in ec_amount:
            ec1 = ec1_list[i]
            ec2 = ec2_list[i]
            strains_c = ec2+(ec1-ec2)*yci
            strains_s = ec2+(ec1-ec2)*ysi
            nc_sum = numpy.sum(concrete(strains_c)*dAc)
            ns_sum = numpy.sum(rebar(strains_s)*omegai)
            mc_sum = numpy.sum(concrete(strains_c)*dAc*(0.5-yci))
            ms_sum = numpy.sum(rebar(strains_s)*omegai*(0.5 - ysi))
            vd = -(nc_sum+ns_sum)
            md = mc_sum+ms_sum
            vd_list.append(vd)
            md_list.append(md)
        data[d1_h]['curves'][omegatot] = {'vd': vd_list, 'md': md_list}
    critical_pairs = []
    for es in [x/10 for x in range(-round(eyd*10),22,2)]:
        critical_pairs.append((ecu2,es))
    for es in [x/10 for x in range(-round(eyd*10),-round(eyd*10),5)]:
        critical_pairs.append((ecu2,es))
    if d1_h==0.05:
        for es in [x/10 for x in range(25,40,5)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(41,76,5)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(84,260,15)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(260,1010,50)]:
            critical_pairs.append((ecu2,-es))
    if d1_h==0.10:
        for es in [x/10 for x in range(25,40,4)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(41,71,6)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(76,140,8)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(140,202,8)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(210,470,20)]:
            critical_pairs.append((ecu2,-es))
    if d1_h==0.15:
        for es in [x/10 for x in range(28,40,4)]:
            critical_pairs.append((ecu2,-es))
        critical_pairs.append((ecu2,-2.5))
        for es in [x/10 for x in range(40,110,5)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(110,160,8)]:
            critical_pairs.append((ecu2,-es))
        for es in [x/10 for x in range(180,280,20)]:
            critical_pairs.append((ecu2,-es))
    if d1_h==0.20:
        critical_pairs.append((ecu2,-2.5))
        critical_pairs.append((ecu2,-2.7))
        critical_pairs.append((ecu2,-3.3))
        critical_pairs.append((ecu2,-4.10))
        critical_pairs.append((ecu2,-4.8))
        critical_pairs.append((ecu2,-5.4))
        critical_pairs.append((ecu2,-7.5))
        critical_pairs.append((ecu2,-6.8))
        critical_pairs.append((ecu2,-6.1))
        critical_pairs.append((ecu2,-8.9))
        critical_pairs.append((ecu2,-8.2))
        critical_pairs.append((ecu2,-11))
        critical_pairs.append((ecu2,-11.8))
        critical_pairs.append((ecu2,-10.3))
        critical_pairs.append((ecu2,-9.70))
        for es in [x/10 for x in range(30,194,7)]:
            critical_pairs.append((ecu2,-es))
    for ec2, es1 in critical_pairs:
        vd_list = []
        md_list = []
        grid = 200
        rsi = 0.5-d1_h
        angles = numpy.linspace(0,2*numpy.pi,bars,endpoint=False)
        ysi = 0.5+rsi*numpy.sin(angles)
        dy = 1/grid
        yci = numpy.linspace(dy/2,1-dy/2,grid)
        bci = 2*numpy.sqrt(0.5**2-(yci-0.5)**2)
        Ac = numpy.pi*(0.5**2)
        dAc = (bci*dy)/Ac
        r = (es1-ec2)/(1-d1_h)
        ec1 = ec2+r
        for omegatot in [0,2.1]:
            omegai = numpy.full(bars,omegatot/bars)
            strains_c = ec2+(ec1-ec2)*yci
            strains_s = ec2+(ec1-ec2)*ysi
            nc_sum = numpy.sum(concrete(strains_c)*dAc)
            ns_sum = numpy.sum(rebar(strains_s)*omegai)
            mc_sum = numpy.sum(concrete(strains_c)*dAc*(0.5-yci))
            ms_sum = numpy.sum(rebar(strains_s)*omegai*(0.5 - ysi))
            vd = -(nc_sum+ns_sum)
            md = mc_sum+ms_sum
            vd_list.append(vd)
            md_list.append(md)
        data[d1_h]['lines'][(ec2, es1)] = {'vd': vd_list, 'md': md_list}
from matplotlib import pyplot
import math
import re
from matplotlib.ticker import MultipleLocator
pyplot.rcParams['svg.fonttype'] = 'path'
pyplot.rcParams['font.family'] = 'Times New Roman'
pyplot.rcParams['font.size'] = 9
pyplot.rcParams['axes.labelsize'] = 9
pyplot.rcParams['xtick.labelsize'] = 9
pyplot.rcParams['ytick.labelsize'] = 9
pyplot.rcParams['legend.fontsize'] = 9
pyplot.rcParams['mathtext.fontset'] = 'cm'
for d1_h in d1_h_list:
    plotdata = data[d1_h]
    pyplot.figure(figsize=(11.2,6.6))
    pyplot.xlim(0,1.2)
    pyplot.ylim(3.4,-4.4)
    ax = pyplot.gca()
    ax = pyplot.gca()
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.figure.canvas.draw()
    for omegatot, vd_md in plotdata['curves'].items():
        omega_int = int(round(omegatot*100))
        is_main_curve = (omega_int%20==0)
        if is_main_curve:
            pyplot.plot(vd_md['md'],vd_md['vd'],linewidth=1,color='k',linestyle='-')
            max_md = max(vd_md['md'])
            max_idx = vd_md['md'].index(max_md)
            vd_nose = vd_md['vd'][max_idx]
            vd_top = vd_md['vd'][-1]
            target_vd = vd_nose+0.50*(vd_top-vd_nose)
            min_diff = float('inf')
            idx = max_idx
            for i in range(max_idx,len(vd_md['vd'])):
                diff = abs(vd_md['vd'][i]-target_vd)
                if diff < min_diff:
                    min_diff = diff
                    idx = i
            x_text, y_text = vd_md['md'][idx], vd_md['vd'][idx]
            idx_prev = max(max_idx,idx-15)
            idx_next = min(len(vd_md['md'])-1,idx+15)
            p1 = ax.transData.transform((vd_md['md'][idx_prev],vd_md['vd'][idx_prev]))
            p2 = ax.transData.transform((vd_md['md'][idx_next],vd_md['vd'][idx_next]))
            dx = p2[0]-p1[0]
            dy = p2[1]-p1[1]
            angle = math.degrees(math.atan2(dy,dx))
            if angle > 90:
                angle -= 180
            elif angle < -90:
                angle += 180
            angle += 3.0
            pyplot.text(x_text,y_text-0.02,'$\omega_{tot}$'+f'={omegatot:.1f}',
                        fontsize=7,
                        rotation=angle,
                        rotation_mode='anchor',
                        verticalalignment='bottom',
                        horizontalalignment='center')
        else:
            pyplot.plot(vd_md['md'],vd_md['vd'],linewidth=0.4,color='silver',linestyle='-',alpha=0.7)
    angle_correction = -2
    for pairs, vd_md in plotdata['lines'].items():
        pyplot.plot(vd_md['md'],vd_md['vd'],linewidth=0.5,color="#000599",dashes=(4, 4))
        nums = re.findall(r"[-+]?(?:\d*\.\d+|\d+)",str(pairs))
        if len(nums) >= 2:
            ec2_val = -float(nums[0])
            es1_val = -float(nums[1])
            text_str = f"$\epsilon_{{c2}}/\epsilon_{{s1}}$={ec2_val:.2f}/{es1_val:.2f}"
        else:
            text_str = f"$\epsilon_{{c2}}/\epsilon_{{s1}}$={pairs}"
        x_start, y_start = vd_md['md'][0], vd_md['vd'][0]
        x_end, y_end = vd_md['md'][-1], vd_md['vd'][-1]
        p1 = ax.transData.transform((x_start,y_start))
        p2 = ax.transData.transform((x_end,y_end))
        dx = p2[0]-p1[0]
        dy = p2[1]-p1[1]
        angle = math.degrees(math.atan2(dy,dx))
        if dx >= 0:
            align_ha = 'left'
            text_str = "  "+text_str
        else:
            align_ha = 'right'
            text_str = text_str+"  " 
        if angle > 90:
            angle -= 180
        elif angle < -90:
            angle += 180
        final_angle = angle+angle_correction
        pyplot.text(x_end,y_end,text_str,
                    fontsize=7,
                    color="#000599",
                    rotation=final_angle,
                    rotation_mode='anchor',
                    verticalalignment='center',
                    horizontalalignment=align_ha)
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    pyplot.grid(True,which='major',linestyle='-',linewidth=0.5,color="#9c9c9c")
    pyplot.savefig(f'plot_{d1_h}.png',format='png',dpi=300,bbox_inches='tight')