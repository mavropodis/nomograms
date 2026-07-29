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
    sxi = 1
    syi = 2
    xs_list = []
    ys_list = []
    for i in range(sxi+1):
            xb = b1_b+i*(1-2*b1_b)/sxi
            for j in range(syi + 1):
                yh = h1_h+j*(1-2*h1_h)/syi
                if i==0 or i==sxi or j==0 or j==syi:
                    xs_list.append(xb)
                    ys_list.append(yh)
    xsi = numpy.array(xs_list)
    ysi = numpy.array(ys_list)
    omegai = omegatot / len(xsi)
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
axial_groups = [
    [0.1, 0.0, -0.1, -0.2],
    [-0.3, -0.4, -0.5, -0.6]
]
omegatot_list = numpy.linspace(0,1.4,71)
pyplot.rcParams['svg.fonttype'] = 'path'
pyplot.rcParams['font.family'] = 'Times New Roman'
pyplot.rcParams['font.size'] = 9
pyplot.rcParams['axes.labelsize'] = 9
pyplot.rcParams['xtick.labelsize'] = 9
pyplot.rcParams['ytick.labelsize'] = 9
pyplot.rcParams['legend.fontsize'] = 9
pyplot.rcParams['mathtext.fontset'] = 'cm'
quadrant_angles = {
    1: 45.0,
    2: 135.0,
    3: -135.0,
    4: -45.0
}
for b1_b in [0.05,0.10,0.15,0.20]:
    h1_h = b1_b
    for group_idx, vexternal_group in enumerate(axial_groups,start=1):
        fig = pyplot.figure(figsize=(11.2,8))
        ax = fig.add_axes([0.05,0.35,0.90,0.636])
        for quadrant, vexternal in enumerate(vexternal_group,start=1):
            vinternal = -vexternal
            target_rad = math.radians(quadrant_angles[quadrant])
            for omegatot in omegatot_list:
                mx, my = biaxial(h1_h,b1_b,omegatot,vinternal)
                if len(mx)>0:
                    mx = numpy.array(mx)
                    my = numpy.array(my)
                    first_octant = mx>=my
                    second_octant = mx<=my
                    if quadrant == 1:
                        m1, m2 = mx, my
                    elif quadrant == 2:
                        m1, m2 = -mx, my
                    elif quadrant == 3:
                        m1, m2 = -mx, -my
                    elif quadrant == 4:
                        m1, m2 = mx, -my
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
                            diff = abs(math.atan2(math.sin(ang-target_rad), math.cos(ang-target_rad)))
                            if diff < min_diff:
                                min_diff = diff
                                best_idx = i
                        x_text, y_text = m1[best_idx], m2[best_idx]
                        idx_prev = max(0, best_idx-2)
                        idx_next = min(len(m1)-1, best_idx+2)
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
                                xy=(x_text,y_text), 
                                xytext=(offset_x,offset_y), 
                                textcoords='offset points',
                                fontsize=7,
                                color='black',
                                rotation=angle,
                                verticalalignment='center',
                                horizontalalignment='center',
                                zorder=3)
            Q = 0.46 
            nu_positions = {
                1: (Q, Q),
                2: (-Q, Q),
                3: (-Q, -Q),
                4: (Q, -Q)
            }
            text_x, text_y = nu_positions[quadrant]
            ax.text(text_x, text_y,
                    f'$\\nu_d = {vexternal}$',
                    horizontalalignment='center',verticalalignment='center',
                    fontsize=8,
                    zorder=5)
        limit = 0.60
        shade_color = "#d8d8d8"
        ax.fill([-limit, 0, 0, -limit], [0, 0, limit, limit], color=shade_color, zorder=0) 
        ax.fill([0, limit, limit, 0], [-limit, -limit, 0, 0], color=shade_color, zorder=0)
        ax.plot([-limit, limit, limit, -limit, -limit], 
                [limit, limit, -limit, -limit, limit], 
                color='black', linewidth=1.2, zorder=4)
        cross_arrow = dict(arrowstyle="<|-|>", color='black', lw=1.5, mutation_scale=12, facecolor='black')
        ax.annotate('', xy=(limit, 0), xytext=(-limit, 0), arrowprops=cross_arrow, annotation_clip=False, zorder=5)
        ax.annotate('', xy=(0, limit), xytext=(0, -limit), arrowprops=cross_arrow, annotation_clip=False, zorder=5)
        off_axis = 0.03
        ax.text(limit + off_axis, 0, '$\\mu_y$', ha='left', va='center', fontsize=11, clip_on=False)
        ax.text(-limit - off_axis, 0, '$\\mu_y$', ha='right', va='center', fontsize=11, clip_on=False)
        ax.text(0, limit + off_axis, '$\\mu_z$', ha='center', va='bottom', fontsize=11, clip_on=False)
        ax.text(0, -limit - off_axis, '$\\mu_z$', ha='center', va='top', fontsize=11, clip_on=False)
        ticks = numpy.round(numpy.arange(-limit, limit + 0.01, 0.1), 1)
        labels = []
        for t in ticks:
            if abs(t) >= limit or abs(t) < 0.05:
                labels.append("")
            else:
                labels.append(f"{abs(t):.1f}")
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_yticklabels(labels, fontsize=9)
        ax.tick_params(top=True, bottom=True, left=True, right=True,
                   labeltop=True, labelbottom=True, labelleft=True, labelright=True,
                   length=0, pad=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.axhline(0, color='black', linewidth=1.5)
        ax.axvline(0, color='black', linewidth=1.5)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.minorticks_on()
        ax.grid(True, which='major', linestyle='-', linewidth=0.3, color="#9c9c9c", zorder=1)
        ax.grid(True, which='minor', linestyle=':', linewidth=0.3, color="#9c9c9c", zorder=1)
        pyplot.savefig(f'plot_{h1_h}_part{group_idx}.png', format='png', dpi=300, bbox_inches='tight')