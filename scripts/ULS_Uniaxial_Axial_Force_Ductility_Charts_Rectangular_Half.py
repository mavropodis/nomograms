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
from matplotlib import ticker
ecu2 = 3.5         #1/1000
Es = 200000        #MPa
fyk = 500          #Mpa
eyk = fyk/Es*1000  #1/1000
gammas = 1.15
eyd = eyk/gammas   #1/1000
eud = 200          #1/1000
def azxi(ec1,ec2):
    if ec1<=0 and ec2<=0:
        a = 0
        xi = 0
    elif ec1<=0 and ec2>0:
        if ec2<=2:
            a = 1/12*ec2*(6-ec2)
        else:
            a = (3*ec2-2)/(3*ec2)
        xi = ec2/(ec2-ec1)
    elif ec1>0 and ec2>0 and ec1!=ec2:
        if ec1<=2 and ec2<=2:
            a = (ec2**2*(6-ec2)-ec1**2*(6-ec1))/(12*(ec2-ec1))
        elif ec1<=2 and ec2>2:
            a = (12*ec2-8-ec1**2*(6-ec1))/(12*(ec2-ec1))
        else:
            a = 1
        xi = 1
    else:
        a = 1
        xi = 1
    return a, xi
def vdinternal(ec1,ec2,omegatot,d1_h):
    r = ec1-ec2
    es1 = ec2+r*(1-d1_h)
    es2 = ec2+r*d1_h
    a, xi = azxi(ec1, ec2)
    s1_fyd = max(-1,min(1,es1/eyd))
    s2_fyd = max(-1,min(1,es2/eyd))
    vd = a*xi+0.5*omegatot*(s1_fyd+s2_fyd)
    return vd
def rury(vd,vd0,vd1,vd2,vd3,omegatot,d1_h):
    error = 1e-5
    iterations = 100
    ru = 0
    if abs(vd-vd0)<1e-6 or abs(vd-vd2)<1e-6:
        ru = 0
    elif vd<vd1:
        rmin = (-eud-ecu2)/(1-d1_h)
        rmax = 0
        for _ in range(iterations):
            rmid = (rmin+rmax)/2
            ec2min = -eud-rmin*(1-d1_h)
            ec1min = -eud+rmin*d1_h
            vdmin = vdinternal(ec1min,ec2min,omegatot,d1_h)-vd
            ec2mid = -eud-rmid*(1-d1_h)
            ec1mid = -eud+rmid*d1_h
            vdmid = vdinternal(ec1mid,ec2mid,omegatot,d1_h)-vd
            if abs(vdmid)<error:
                break
            if vdmin*vdmid>0:
                rmin = rmid
            else:
                rmax = rmid
        ru = rmid
    elif vd<vd2:
        rmin = (-eud-ecu2)/(1-d1_h)
        rmax = 0
        for _ in range(iterations):
            rmid = (rmin+rmax)/2
            ec2min = ecu2
            ec1min = ec2min+rmin
            vdmin = vdinternal(ec1min,ec2min,omegatot,d1_h)-vd
            ec2mid = ecu2
            ec1mid = ec2mid+rmid
            vdmid = vdinternal(ec1mid,ec2mid,omegatot,d1_h)-vd
            if abs(vdmid)<error:
                break
            if vdmin*vdmid>0:
                rmin = rmid
            else:
                rmax = rmid
        ru = rmid
    ry = ru
    index = 0
    if vd<=vd3:
        rmin = (-eyd-ecu2)/(1-d1_h)
        rmax = 0
        for _ in range(iterations):
            rmid = (rmin+rmax)/2
            ec2min = -eyd-rmin*(1-d1_h)
            ec1min = -eyd+rmin*d1_h
            vdmin = vdinternal(ec1min,ec2min,omegatot,d1_h)-vd
            ec2mid = -eyd-rmid*(1-d1_h)
            ec1mid = -eyd+rmid*d1_h
            vdmid = vdinternal(ec1mid,ec2mid,omegatot,d1_h)-vd
            if abs(vdmid)<error:
                break
            if vdmin*vdmid>0:
                rmin = rmid
            else:
                rmax = rmid
        ry = rmid
        index = 1
    return ru, ry, index
omegatot_list = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
d1_h_list = [0.05,0.10,0.15,0.20]
pyplot.rcParams['font.family'] = 'Times New Roman'
pyplot.rcParams['font.size'] = 9
pyplot.rcParams['mathtext.fontset'] = 'cm'
fig = pyplot.figure(figsize=(9.5,6.1))
distinct_colors = [
"#7a1400",
"#c60000",
"#000AC5",
"#e76c00",
"#c1c100",
"#00bb00",
"#0073c6",
"#6e00d5",
"#bc0084",
"#00b3c7",
"#e8c500"
]
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], wspace=0, hspace=0)
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[0, 1], sharex=ax0, sharey=ax0)
ax2 = fig.add_subplot(gs[1, 0], sharex=ax0, sharey=ax0)
ax3 = fig.add_subplot(gs[1, 1], sharex=ax0, sharey=ax0)
axes = [ax0, ax1, ax2, ax3]
lines_for_legend = []
labels_for_legend = []
for idx, (ax, d1_h) in enumerate(zip(axes, d1_h_list)):
    for c_idx, omegatot in enumerate(omegatot_list):
        vd0 = -omegatot
        r1 = (-eud-ecu2)/(1-d1_h)
        ec21 = ecu2
        ec11 = ec21+r1
        vd1 = vdinternal(ec11,ec21,omegatot,d1_h)
        ec12 = ec22 = 3.5
        vd2 = vdinternal(ec12,ec22,omegatot,d1_h)
        r3 = (-eyd-ecu2)/(1-d1_h)
        ec23 = ecu2
        ec13 = ec23+r3
        vd3 = vdinternal(ec13,ec23,omegatot,d1_h)
        vd_array = numpy.linspace(-0.2,0.5,500)
        mphi_data = []
        vd_data = []
        for vd in vd_array:
            ru, ry, index = rury(vd,vd0,vd1,vd2,vd3,omegatot,d1_h)
            if index==1 and abs(ry)>1e-8 and abs(ru)>1e-8:
                mphi = abs(ru)/abs(ry)
            else:
                mphi = 1
            mphi_data.append(mphi)
            vd_data.append(-vd)
        line, = ax.plot(mphi_data, vd_data, color=distinct_colors[c_idx], linewidth=1.5)
        if idx == 0:
            lines_for_legend.append(line)
            labels_for_legend.append(f'$\\omega_{{tot}} = {omegatot}$')
    ax.set_xscale('log')
    ax.set_xticks([1,2,5,10,20,50,100,200])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=numpy.arange(2, 10.0)))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.grid(True, which='major', linestyle='-', linewidth=0.7, color='gray')
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, color='darkgray', alpha=0.7)
    if idx == 3:
        ax.legend(lines_for_legend, labels_for_legend, loc='upper right', bbox_to_anchor=(0.98, 0.98), fontsize=9)
ax0.tick_params(labelbottom=False)
ax1.tick_params(labelbottom=False, labelleft=False)
ax3.tick_params(labelleft=False)
ax0.set_xlim(left=1, right=250) 
ax0.set_ylim(bottom=0.18, top=-0.5)
pyplot.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.08)
pyplot.savefig(f'plot.png', format='png', dpi=300, bbox_inches='tight')