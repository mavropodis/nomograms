# ==============================================================================
# Copyright (c) 2026 I. Mavropodis, Ch. Zeris. All rights reserved.
# 
# This source code is part of the Diploma Thesis developed at the 
# National Technical University of Athens (NTUA).
# It is provided for review and academic purposes only. 
# Unauthorized copying, modification, or commercial use is strictly prohibited.
# ==============================================================================
import numpy
ecu2 = 3.5         #1/1000
Es = 200000        #MPa
fyk = 500          #Mpa
eyk = fyk/Es*1000  #1/1000
gammas = 1.15
eyd = eyk/gammas   #1/1000
eud = 200          #1/1000
def az(ec1,ec2):
    if ec1<=0 and ec2<=0:
        a = 0
        z = 0
    elif ec1<=0 and ec2>0:
        if ec2<=2:
            a = 1/12*ec2*(6-ec2)
            z = (8-ec2)/(4*(6-ec2))
        else:
            a = (3*ec2-2)/(3*ec2)
            z = (ec2*(3*ec2-4)+2)/(2*ec2*(3*ec2-2))
    elif ec1>0 and ec2>0 and ec1!=ec2:
        if ec1<=2 and ec2<=2:
            a = (ec2**2*(6-ec2)-ec1**2*(6-ec1))/(12*(ec2-ec1))
            z = (8*ec2**3-ec2**4-24*ec1**2*ec2+4*ec1**3*ec2-3*ec1**4+16*ec1**3)/(4*(6*ec2**2-ec2**3-6*ec1**2+ec1**3)*(ec2-ec1))
        elif ec1<=2 and ec2>2:
            a = (12*ec2-8-ec1**2*(6-ec1))/(12*(ec2-ec1))
            z = (4*ec1**3*ec2-24*ec1**2*ec2-3*ec1**4+16*ec1**3+24*ec2**2-32*ec2+16)/(4*(ec1**3-6*ec1**2+12*ec2-8)*(ec2-ec1))
        else:
            a = 1
            z = 0.5
    else:
        a = 1
        z = 0.5
    return a, z
def vminternal(Phi,xi,omegatot,d1_h):
    ec1 = Phi*(xi-1)
    ec2 = Phi*xi
    a, z = az(ec1, ec2)
    if ec1>0:
        xic = 1
    else:
        xic = xi
    n = 9
    top_bottom = 5/24
    middle = 2/24
    ns_sum = 0
    ms_sum = 0
    for i in range(1,n+1):
        if i==1 or i==n:
            omegai = top_bottom*omegatot
        else:
            omegai = middle*omegatot
        di = d1_h+(i-1)/(n-1)*(1-2*d1_h)
        esi = Phi*(xi-di)
        si_fyd = max(-1,min(1,esi/eyd))
        ns_sum += si_fyd*omegai
        ms_sum += si_fyd*omegai*(0.5-di)
    vinternal = a*xic+ns_sum
    minternal = a*xic*(0.5-z*xic)+ms_sum
    return vinternal, minternal
def ximinternalmid(Phi,vexternal,omegatot,d1_h):
    ximin = -2
    ximax = 10
    error = 1e-5
    iterations = 100
    vinternalmin, minternalmin = vminternal(Phi, ximin, omegatot, d1_h)
    vinternalmax, minternalmax = vminternal(Phi, ximax, omegatot, d1_h)
    dvmin = vinternalmin-vexternal
    dvmax = vinternalmax-vexternal
    if dvmin*dvmax>0: 
        return None, None
    for _ in range(iterations):
        ximid = (ximin+ximax)/2
        vinternalmid, minternalmid = vminternal(Phi, ximid, omegatot, d1_h)
        dvmid = vinternalmid-vexternal
        if abs(dvmid)<error:
            return ximid, minternalmid
        if dvmid*dvmin<0:
            ximax = ximid
            dvmax = dvmid
        else:
            ximin = ximid
            dvmin = dvmid
    return ximid, minternalmid
import numpy
from matplotlib import pyplot
import matplotlib.ticker as ticker
d1_h_list = [0.05, 0.10, 0.15, 0.20]
omegatot_pages = [
    [0.10, 0.15, 0.20, 0.25],
    [0.30, 0.35, 0.40, 0.45],
    [0.50, 0.55, 0.60, 0.65]
]
vexternal_list = [-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40]
Phi_list = numpy.linspace(0.01, 50, 700)
distinct_colors = [
"#87007383",
"#0087C5",
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
pyplot.rcParams['font.family'] = 'Times New Roman'
pyplot.rcParams['font.size'] = 11
pyplot.rcParams['mathtext.fontset'] = 'cm'
box_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.8, alpha=0.9)
for d1_h in d1_h_list:
    for page_idx, omegatot_list in enumerate(omegatot_pages):
        fig = pyplot.figure(figsize=(11.69, 8.27))
        gs = fig.add_gridspec(2, 2, wspace=0.0, hspace=0.0)
        axes = []
        for i in range(2):
            for j in range(2):
                if not axes:
                    ax = fig.add_subplot(gs[i, j])
                else:
                    ax = fig.add_subplot(gs[i, j], sharex=axes[0], sharey=axes[0])
                axes.append(ax)
        page_max_phi = 0.0
        page_max_md = 0.0
        lines_for_legend = []
        labels_for_legend = []
        for idx, (ax, omegatot) in enumerate(zip(axes, omegatot_list)):
            for c_idx, vexternal in enumerate(vexternal_list):
                Phi_data = []
                md_data = []
                yield_point = None
                for Phi in Phi_list:
                    xi, md = ximinternalmid(Phi, vexternal, omegatot, d1_h)
                    if xi is not None:
                        ec1 = Phi * (xi - 1)
                        ec2 = Phi * xi
                        es1 = Phi * (xi + d1_h - 1)
                        es2 = Phi * (xi - d1_h)
                        if ec2 > ecu2 or (abs(es1) > eud and es1 < 0) or (abs(es2) > eud and es2 < 0):
                            break
                        Phi_data.append(Phi)
                        md_data.append(md)
                        if yield_point is None and es1 <= -eyd:
                            yield_point = (Phi, md)
                if Phi_data:
                    Phi_data.insert(0, 0)
                    md_data.insert(0, 0)
                    color_to_use = distinct_colors[c_idx % len(distinct_colors)]
                    line, = ax.plot(Phi_data, md_data, color=color_to_use, linewidth=1.5)
                    last_phi = Phi_data[-1]
                    last_md = md_data[-1]
                    page_max_phi = max(page_max_phi, last_phi)
                    page_max_md = max(page_max_md, max(md_data))
                    if idx == 0:
                        lines_for_legend.append(line)
                        labels_for_legend.append(f'$\\nu_d = {vexternal:.2f}$')
                    if yield_point:
                        ax.plot(yield_point[0], yield_point[1], marker='o', markersize=1.5, 
                                color=color_to_use, markeredgewidth=1)
            ax.grid(True, which='major', linestyle='-', linewidth=0.7, color='gray')
            ax.grid(True, which='minor', linestyle=':', linewidth=0.5, color='lightgray')
            ax.minorticks_on()
            ax.text(0.95, 0.95, f'$\\omega_{{tot}} = {omegatot:.2f}$', transform=ax.transAxes, 
                    fontsize=12, fontweight='bold', ha='right', va='top', bbox=box_props)
            if idx < 2:
                ax.tick_params(labelbottom=False)
            if idx % 2 != 0:
                ax.tick_params(labelleft=False)
            if idx == 3:
                ax.legend(lines_for_legend, labels_for_legend, loc='lower right', 
                          fontsize=10, framealpha=1.0, edgecolor='black')
        if page_max_phi > 0 and page_max_md > 0:
            axes[0].set_xlim(left=0, right=page_max_phi * 1.05)
            axes[0].set_ylim(bottom=0, top=page_max_md * 1.10)
        for ax in axes:
            ax.xaxis.set_major_locator(ticker.MaxNLocator(prune='upper'))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(prune='upper'))
        if page_max_phi > 0 and page_max_md > 0:
            axes[0].set_xlim(left=0, right=page_max_phi * 1.15)
            axes[0].set_ylim(bottom=0, top=page_max_md * 1.10)
        pyplot.subplots_adjust(left=0.05, right=0.98, top=0.95, bottom=0.08)
        pyplot.savefig(f'plot_d1h_{d1_h}_page_{page_idx + 1}.svg', format='svg', bbox_inches='tight')
        pyplot.close(fig)