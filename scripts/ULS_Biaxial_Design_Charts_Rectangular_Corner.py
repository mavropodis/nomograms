# ==============================================================================
# Copyright (c) 2026 I. Mavropodis, Ch. Zeris. All rights reserved.
# 
# This source code is part of the Diploma Thesis developed at the 
# National Technical University of Athens (NTUA).
# It is provided for review and academic purposes only. 
# Unauthorized copying, modification, or commercial use is strictly prohibited.
# ==============================================================================
import numpy
from matplotlib import pyplot
ec2, ecu2 = 2.0, 3.5         #1/1000
Es = 200000                  #MPa
fyk = 500                    #Mpa
eyk = fyk/Es*1000            #1/1000
gammas = 1.15
eyd = eyk/gammas             #1/1000
eud = 200                    #1/1000
def concrete(strains):
    sc_fcd = numpy.zeros_like(strains)
    parabola = (strains>0) & (strains<=ec2)
    sc_fcd[parabola] = 1-(1-strains[parabola]/ec2)**2
    rectangular = (strains>ec2)
    sc_fcd[rectangular] = 1
    return sc_fcd
def rebar(strains):
    ss_fyd = numpy.zeros_like(strains)
    elastic = numpy.abs(strains)<=eyd
    ss_fyd[elastic] = strains[elastic]/eyd
    plastic = (numpy.abs(strains)>eyd)
    ss_fyd[plastic] = numpy.sign(strains[plastic])
    return ss_fyd
def biaxial(h1_h, b1_b, omegatot, vinternal):
    xs1, ys1 = b1_b, h1_h
    xs2, ys2 = 1-b1_b, h1_h
    xs3, ys3 = 1-b1_b, 1-h1_h
    xs4, ys4 = b1_b, 1-h1_h
    xsi = numpy.array([xs1, xs2, xs3, xs4])
    ysi = numpy.array([ys1, ys2, ys3, ys4])
    omegai = omegatot/4
    grid = 80
    step = 1/grid
    coordinates = numpy.linspace(step/2, 1-step/2, grid)
    xci, yci = numpy.meshgrid(coordinates,coordinates)
    xci, yci = xci.flatten(), yci.flatten()
    dAc = step*step
    mx_list = []
    my_list = []
    thetas = numpy.linspace(0.00001,89.9999,500)
    for theta in thetas:
        theta = numpy.radians(theta)
        yci_tonos = xci*numpy.cos(theta)+yci*numpy.sin(theta)
        yci_tonos_max = numpy.max(yci_tonos)
        ysi_tonos = xsi*numpy.cos(theta)+ysi*numpy.sin(theta)
        found = False
        rmin = -10000
        rmax = 5000000
        error = 1e-5
        iterations = 700
        for _ in range(iterations):
            rmid = (rmin+rmax)/2
            e0max = ecu2-rmid*yci_tonos_max
            eci = e0max+rmid*yci_tonos
            esi = e0max+rmid*ysi_tonos
            vc = numpy.sum(concrete(eci)*dAc)
            vs = numpy.sum(rebar(esi)*omegai)
            vinternalmid = vc+vs
            dvmid = vinternalmid-vinternal
            if abs(dvmid)<error:
                mcx = numpy.sum(concrete(eci)*dAc*(yci-0.50))
                mcy = numpy.sum(concrete(eci)*dAc*(xci-0.50))
                msx = numpy.sum(rebar(esi)*omegai*(ysi-0.50))
                msy = numpy.sum(rebar(esi)*omegai*(xsi-0.50))
                mx_list.append(mcx+msx)
                my_list.append(mcy+msy)
                found = True
                break
            if dvmid>0:
                rmin = rmid
            else:
                rmax = rmid
            if not found:
                pass
    return mx_list, my_list
import math
vexternal_list = [0.1,0.0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6]
omegatot_list = numpy.linspace(0,1.4,71)
pyplot.rcParams['svg.fonttype'] = 'path'
pyplot.rcParams['font.family'] = 'Times New Roman'
pyplot.rcParams['font.size'] = 9
pyplot.rcParams['axes.labelsize'] = 9
pyplot.rcParams['xtick.labelsize'] = 9
pyplot.rcParams['ytick.labelsize'] = 9
pyplot.rcParams['legend.fontsize'] = 9
pyplot.rcParams['mathtext.fontset'] = 'cm'
octant_angles = {
    1: 67.5,
    2: 22.5,
    3: -22.5,
    4: -67.5,
    5: -112.5,
    6: -157.5,
    7: 157.5,
    8: 112.5
}
for b1_b in [0.05,0.10,0.15,0.20]:
    fig = pyplot.figure(figsize=(11.2,8))
    ax = fig.add_axes([0.05,0.35,0.90,0.636])
    h1_h = b1_b
    for octant, vexternal in enumerate(vexternal_list, start=1):
        vinternal = -vexternal
        target_rad = math.radians(octant_angles[octant])
        for omegatot in omegatot_list:
            mx, my = biaxial(h1_h,b1_b,omegatot,vinternal)
            if len(mx)>0:
                mx = numpy.array(mx)
                my = numpy.array(my)
                first_octant = mx>=my
                second_octant = mx<=my
                if octant==1: m1, m2 = mx[second_octant], my[second_octant]
                elif octant==2: m1, m2 = mx[first_octant], my[first_octant]
                elif octant==3: m1, m2 = mx[first_octant], -my[first_octant]
                elif octant==4: m1, m2 = mx[second_octant], -my[second_octant]
                elif octant==5: m1, m2 = -mx[second_octant], -my[second_octant]
                elif octant==6: m1, m2 = -mx[first_octant], -my[first_octant]
                elif octant==7: m1, m2 = -mx[first_octant], my[first_octant]
                elif octant==8: m1, m2 = -mx[second_octant], my[second_octant]
                major_omegatot = round(omegatot*100)%20==0
                linewidth = 1.0 if major_omegatot else 0.5
                color = 'black' if major_omegatot else "#7D7D7D"
                alpha = 1.0 if major_omegatot else 0.7
                ax.plot(m1,m2,color=color,linewidth=linewidth,alpha=alpha,zorder=2)
                if major_omegatot and len(m1) > 1:
                    best_idx = 0
                    min_diff = float('inf')
                    for i in range(len(m1)):
                        ang = math.atan2(m2[i], m1[i])
                        diff = abs(math.atan2(math.sin(ang-target_rad),math.cos(ang-target_rad)))
                        if diff < min_diff:
                            min_diff = diff
                            best_idx = i
                    x_text, y_text = m1[best_idx], m2[best_idx]
                    idx_prev = max(0,best_idx-2)
                    idx_next = min(len(m1)-1,best_idx+2)
                    if idx_prev != idx_next:
                        dx = m1[idx_next]-m1[idx_prev]
                        dy = m2[idx_next]-m2[idx_prev]
                        angle = math.degrees(math.atan2(dy,dx))
                    else:
                        angle = 0
                    if angle > 90:   
                        angle -= 180
                    elif angle < -90: 
                        angle += 180
                    shift_points = 6
                    offset_x = math.cos(target_rad)*shift_points
                    offset_y = math.sin(target_rad)*shift_points
                    if omegatot!=1.0:
                        text_omegatot = f'{omegatot:.1f}'
                    else:
                        text_omegatot = '$\\omega_{tot}$='+f'{omegatot:.1f}'
                    ax.annotate(text_omegatot, 
                            xy=(x_text, y_text), 
                            xytext=(offset_x, offset_y), 
                            textcoords='offset points',
                            fontsize=7,
                            color='black',
                            rotation=angle,
                            verticalalignment='center',
                            horizontalalignment='center',
                            zorder=3)
        L = 0.52
        S = 0.42
        nu_positions = {
            1: (S, L),
            2: (L, S),
            3: (L, -S),
            4: (S, -L),
            5: (-S, -L),
            6: (-L, -S),
            7: (-L, S),
            8: (-S, L)}
        text_x, text_y = nu_positions[octant]
        ax.text(text_x, text_y,
                f'$\\nu_d = {vexternal}$',
                horizontalalignment='center',verticalalignment='center',
                fontsize=8,
                zorder=5)
    limit = 0.60
    shade_color = "#d8d8d8"
    ax.fill([0, limit, 0], [0, limit, limit], color=shade_color, zorder=0)
    ax.fill([0, limit, limit], [0, 0, -limit], color=shade_color, zorder=0)
    ax.fill([0, -limit, 0], [0, -limit, -limit], color=shade_color, zorder=0)
    ax.fill([0, -limit, -limit], [0, 0, limit], color=shade_color, zorder=0)
    arrow_style = dict(arrowstyle="<|-|>", color='black', lw=1.2, mutation_scale=10, facecolor='black')
    thick_arrow = dict(arrowstyle="<|-|>", color='black', lw=1.5, mutation_scale=10, facecolor='black')
    ax.plot([-limit,limit], [-limit,limit], color='black', linestyle='-', linewidth=1.5, zorder=4) 
    ax.plot([-limit,limit], [limit,-limit], color='black', linestyle='-', linewidth=1.5, zorder=4)
    ax.annotate('', xy=(limit,0), xytext=(-limit,0), arrowprops=thick_arrow, annotation_clip=False, zorder=4)
    ax.annotate('', xy=(0,limit), xytext=(0,-limit), arrowprops=thick_arrow, annotation_clip=False, zorder=4)
    ax.annotate('', xy=(limit,limit), xytext=(-limit,limit), arrowprops=arrow_style, annotation_clip=False, zorder=4)
    ax.annotate('', xy=(limit,-limit), xytext=(-limit,-limit), arrowprops=arrow_style, annotation_clip=False, zorder=4)
    ax.annotate('', xy=(limit,limit), xytext=(limit,-limit), arrowprops=arrow_style, annotation_clip=False, zorder=4)
    ax.annotate('', xy=(-limit,limit), xytext=(-limit,-limit), arrowprops=arrow_style, annotation_clip=False, zorder=4)
    off1 = 0.025
    ax.text(limit+off1, 0, '$\\mu_1$', ha='left', va='center', fontsize=10, clip_on=False)
    ax.text(-limit-off1, 0, '$\\mu_1$', ha='right', va='center', fontsize=10, clip_on=False)
    ax.text(0, limit+off1, '$\\mu_1$', ha='center', va='bottom', fontsize=10, clip_on=False)
    ax.text(0, -limit-off1, '$\\mu_1$', ha='center', va='top', fontsize=10, clip_on=False)
    off2 = 0.025
    ax.text(limit+off2, limit, '$\\mu_2$', ha='left', va='center', fontsize=10, clip_on=False)
    ax.text(limit, limit+off2, '$\\mu_2$', ha='center', va='bottom', fontsize=10, clip_on=False)
    ax.text(-limit-off2, limit, '$\\mu_2$', ha='right', va='center', fontsize=10, clip_on=False)
    ax.text(-limit, limit+off2, '$\\mu_2$', ha='center', va='bottom', fontsize=10, clip_on=False)
    ax.text(limit+off2, -limit, '$\\mu_2$', ha='left', va='center', fontsize=10, clip_on=False)
    ax.text(limit, -limit-off2, '$\\mu_2$', ha='center', va='top', fontsize=10, clip_on=False)
    ax.text(-limit-off2, -limit, '$\\mu_2$', ha='right', va='center', fontsize=10, clip_on=False)
    ax.text(-limit, -limit-off2, '$\\mu_2$', ha='center', va='top', fontsize=10, clip_on=False)
    ticks = numpy.round(numpy.arange(-limit, limit+0.01, 0.1), 1)
    labels = []
    for t in ticks:
        if abs(t) >= limit or abs(t) < 0.05:
            labels.append("")
        else:
            labels.append(f"{abs(t):.1f}")
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(labels,fontsize=9)
    ax.set_yticklabels(labels,fontsize=9)
    ax.tick_params(top=True,bottom=True,left=True,right=True,
               labeltop=True,labelbottom=True,labelleft=True,labelright=True,
               length=0,pad=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.axhline(0,color='black',linewidth=1.5)
    ax.axvline(0,color='black',linewidth=1.5)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xlim(-limit,limit)
    ax.set_ylim(-limit,limit)
    ax.minorticks_on()
    ax.grid(True,which='major',linestyle='-',linewidth=0.3,color="#9c9c9c",zorder=1)
    ax.grid(True,which='minor',linestyle=':',linewidth=0.3,color="#9c9c9c",zorder=1)
    pyplot.savefig(f'plot_{h1_h}.png',format='png',dpi=300,bbox_inches='tight')