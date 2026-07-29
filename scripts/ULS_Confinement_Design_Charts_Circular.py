# ==============================================================================
# Copyright (c) 2026 I. Mavropodis, Ch. Zeris. All rights reserved.
# 
# This source code is part of the Diploma Thesis developed at the 
# National Technical University of Athens (NTUA).
# It is provided for review and academic purposes only. 
# Unauthorized copying, modification, or commercial use is strictly prohibited.
# ==============================================================================
from matplotlib import pyplot
import math
from matplotlib.ticker import MultipleLocator
pyplot.rcParams['font.family'] = 'Times New Roman'
pyplot.rcParams['font.size'] = 9
pyplot.rcParams['axes.labelsize'] = 9
pyplot.rcParams['xtick.labelsize'] = 9
pyplot.rcParams['ytick.labelsize'] = 9
pyplot.rcParams['legend.fontsize'] = 9
pyplot.rcParams['mathtext.fontset'] = 'cm'
import numpy
eps2, ecu = 2.0, 3.5         #1/1000
Es = 200000                  #MPa
fyk = 500                    #Mpa
eyk = fyk/Es*1000            #1/1000
gammas = 1.15
eyd = eyk/gammas             #1/1000
eud = 200                    #1/1000
bars = 8
d1_h_list = [0.02,0.04,0.08,0.10]
omegaconf_list = numpy.arange(0,1.25,0.05)
omegatot_list = [0.2,0.4,0.6,0.8,1.0]
kconf_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
kaggregate = 1.0
def concrete(strains, ec2c, ecuc, fc):
    sc_fcd = numpy.zeros_like(strains)
    parabola = (strains > 0) & (strains <= ec2c)
    sc_fcd[parabola] = fc * (1 - (1 - strains[parabola]/ec2c)**2)
    rectangular = (strains > ec2c) & (strains <= ecuc)
    sc_fcd[rectangular] = fc
    return sc_fcd
def rebar(strains):
    ss_fyd = numpy.zeros_like(strains)
    elastic = numpy.abs(strains)<=eyd
    ss_fyd[elastic] = strains[elastic]/eyd
    plastic = (numpy.abs(strains)>eyd)
    ss_fyd[plastic] = numpy.sign(strains[plastic])
    return ss_fyd
for d1_h in d1_h_list:
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
        data = {}
        omegai = numpy.full(bars, omegatot / bars)
        for kconf in kconf_list:
            data[kconf]={'curves':{}}
            for omegaconf in omegaconf_list:
                ec1_list = []
                es1_list = []
                es2_list = []
                ec2_list = []
                if omegaconf<=0.6:
                    Df = 4*omegaconf*kaggregate
                else:
                    Df = 3.5*omegaconf**0.75*kaggregate
                fc = 1+kconf*Df
                ec2c = eps2*(1+5*Df)
                ecuc = ecu+200*omegaconf
                step = 0.001
                multiplier = int(1/step)
                r_start_1 = 0
                r_end_1 = (-eud-ecuc)/(1-d1_h)
                start_1 = int(round(r_start_1*multiplier))
                finish_1 = int(round(r_end_1*multiplier))
                for r in range(start_1,finish_1-1,-1):
                    r = r/float(multiplier)
                    es1i = -eud
                    ec1i = es1i+r*d1_h
                    ec2i = es1i-r*(1-d1_h)
                    ec1_list.append(ec1i)
                    ec2_list.append(ec2i)
                r_start_2 = (-eud-ecuc)/(1-d1_h)
                r_end_2 = 0
                start_2 = int(round(r_start_2*multiplier))
                finish_2 = int(round(r_end_2*multiplier))
                for r in range(start_2,finish_2+1,1):
                    r = r/float(multiplier)
                    ec2i = ecuc
                    ec1i = ec2i+r
                    ec1_list.append(ec1i)
                    ec2_list.append(ec2i)
                vd_list = []
                md_list = []
                ec_amount = range(len(ec2_list))
                for i in ec_amount:
                    ec1 = ec1_list[i]
                    ec2 = ec2_list[i]
                    strains_c = ec2+(ec1-ec2)*yci
                    strains_s = ec2+(ec1-ec2)*ysi
                    nc_sum = numpy.sum(concrete(strains_c, ec2c, ecuc, fc)*dAc)
                    ns_sum = numpy.sum(rebar(strains_s)*omegai)
                    mc_sum = numpy.sum(concrete(strains_c, ec2c, ecuc, fc)*dAc*(0.5-yci))
                    ms_sum = numpy.sum(rebar(strains_s)*omegai*(0.5 - ysi))
                    vd = -(nc_sum+ns_sum)
                    md = mc_sum+ms_sum
                    vd_list.append(vd)
                    md_list.append(md)
                data[kconf]['curves'][omegaconf] = {'vd': vd_list, 'md': md_list}
        for page in range(2):
            fig, axes = pyplot.subplots(nrows=2, ncols=3, figsize=(16, 10), sharex=True, sharey=True)
            fig.subplots_adjust(wspace=0.0, hspace=0.0)
            axs = axes.flatten()
            for i in range(6):
                ax = axs[i]
                if i == 5:
                    ax.axis('off')
                    continue
                idx = page * 5 + i
                k_conf_current = kconf_list[idx]
                plotdata = data[k_conf_current]
                ax.set_xlim(0, 1.0)
                ax.set_ylim(0, -2.0)
                ax.xaxis.set_major_locator(MultipleLocator(0.1))
                ax.yaxis.set_major_locator(MultipleLocator(0.2))
                for omega_conf, vd_md in plotdata['curves'].items():
                    if len(vd_md['md']) == 0: continue
                    omega_int = int(round(omega_conf * 100))
                    is_main_curve = (omega_int % 10 == 0)
                    if is_main_curve:
                        line_color = 'k'
                        line_width = 1.0 
                        ax.plot(vd_md['md'], vd_md['vd'], linewidth=line_width, color=line_color, linestyle='-')
                    else:
                        ax.plot(vd_md['md'], vd_md['vd'], linewidth=0.5, color='grey', linestyle='-', alpha=0.8)
                        continue 
                    fig.canvas.draw()
                    max_md = max(vd_md['md'])
                    max_idx = vd_md['md'].index(max_md)
                    vd_nose = vd_md['vd'][max_idx]
                    top_idx = max_idx
                    while top_idx < len(vd_md['md']) - 1:
                        if vd_md['md'][top_idx + 1] < 1e-3:
                            break
                        top_idx += 1
                    vd_top = vd_md['vd'][top_idx]
                    label_step = omega_int // 10
                    is_alternate = (label_step % 2 != 0)
                    placement_ratio = 0.60 if is_alternate else 0.40
                    target_vd = vd_nose + placement_ratio * (vd_top - vd_nose)
                    min_diff = float('inf')
                    idx_text = max_idx
                    for j in range(max_idx, top_idx + 1):
                        diff = abs(vd_md['vd'][j] - target_vd)
                        if diff < min_diff:
                            min_diff = diff
                            idx_text = j
                    x_text, y_text = vd_md['md'][idx_text], vd_md['vd'][idx_text]
                    if y_text < -1.95:
                        placement_ratio = 0.25
                        target_vd = vd_nose + placement_ratio * (vd_top - vd_nose)
                        min_diff = float('inf')
                        for j in range(max_idx, top_idx + 1):
                            diff = abs(vd_md['vd'][j] - target_vd)
                            if diff < min_diff:
                                min_diff = diff
                                idx_text = j
                        x_text, y_text = vd_md['md'][idx_text], vd_md['vd'][idx_text]
                        if y_text < -1.95:
                            continue
                    idx_prev = max(max_idx, idx_text - 15)
                    idx_next = min(top_idx, idx_text + 15)
                    if idx_prev >= idx_next:
                        continue
                    p1 = ax.transData.transform((vd_md['md'][idx_prev], vd_md['vd'][idx_prev]))
                    p2 = ax.transData.transform((vd_md['md'][idx_next], vd_md['vd'][idx_next]))
                    dx = p2[0] - p1[0]
                    dy = p2[1] - p1[1]
                    if dx == 0 and dy == 0:
                        continue
                    raw_angle = math.degrees(math.atan2(dy, dx))
                    if abs(raw_angle) < 15 or abs(raw_angle - 180) < 15 or abs(raw_angle + 180) < 15:
                        continue
                    angle = raw_angle
                    if angle > 90: angle -= 180
                    elif angle < -90: angle += 180
                    angle += 3.0
                    label_text = f'$\\omega_{{conf}}$={omega_conf:.1f}'
                    if omega_conf == 0.0: label_text = '$\\omega_{{conf}}=0.0$'
                    ax.text(x_text, y_text - 0.02, label_text,
                            fontsize=7, rotation=angle, rotation_mode='anchor',
                            verticalalignment='bottom', horizontalalignment='center', color='k', clip_on=True)
                ax.grid(True, which='major', linestyle='-', linewidth=0.5, color="#9c9c9c")
            for i in range(6):
                if i == 5: continue
                ax = axs[i]
                if i % 3 != 0:
                    ax.tick_params(axis='y', labelleft=False)
                else:
                    ax.tick_params(axis='y', labelleft=True)
                if i in [2, 3, 4]:
                    ax.tick_params(axis='x', labelbottom=True)
                else:
                    ax.tick_params(axis='x', labelbottom=False)
                if i % 3 != 0:
                    for tick in ax.xaxis.get_major_ticks():
                        if math.isclose(tick.get_loc(), 0.0, abs_tol=1e-3):
                            tick.label1.set_visible(False)
                if i == 0:
                    for tick in ax.yaxis.get_major_ticks():
                        if math.isclose(tick.get_loc(), 0.0, abs_tol=1e-3):
                            tick.label1.set_visible(False)
            filename = f'Assessment_omegatot_{omegatot:.2f}_d1h_{d1_h:.2f}_Page_{page+1}.png'
            pyplot.savefig(filename, format='png', dpi=300, bbox_inches='tight')
            pyplot.close(fig)