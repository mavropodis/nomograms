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
eps2 , ecu2 = 2.0, 3.5        #1/1000
Es = 200000                   #MPa
fyk = 500                     #MPa
eyk = fyk/Es*1000             #1/1000
gammas = 1.15
eyd = eyk/gammas              #1/1000
eud = 200                     #1/1000
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
def vminternal(ec1, ec2, omegatot, d1_h):
    strains_c = ec2 + (ec1 - ec2) * yci
    sc = concrete(strains_c)
    nc_sum = numpy.sum(sc * dAc)
    mc_sum = numpy.sum(sc * dAc * (0.5 - yci))
    strains_s = ec2 + (ec1 - ec2) * ysi
    ss = rebar(strains_s)
    omegai = omegatot / bars
    ns_sum = numpy.sum(ss * omegai)
    ms_sum = numpy.sum(ss * omegai * (0.5 - ysi))
    vinternal = nc_sum + ns_sum
    minternal = mc_sum + ms_sum
    return vinternal, minternal
def ximinternalmid(Phi,vexternal,omegatot,d1_h):
    ximin = -2
    ximax = 2
    error = 1e-5
    iterations = 100
    ec2_min = Phi*ximin
    ec1_min = Phi*(ximin-1)
    vinternalmin, minternalmin = vminternal(ec1_min, ec2_min, omegatot, d1_h)
    ec2_max = Phi*ximax
    ec1_max = Phi*(ximax-1)
    vinternalmax, minternalmax = vminternal(ec1_max, ec2_max, omegatot, d1_h)
    dvmin = vinternalmin-vexternal
    dvmax = vinternalmax-vexternal
    if dvmin*dvmax>0:
        return None, None
    for _ in range(iterations):
        ximid = (ximin+ximax)/2
        ec2_mid = Phi*ximid
        ec1_mid = Phi*(ximid-1)
        vinternalmid, minternalmid = vminternal(ec1_mid, ec2_mid, omegatot, d1_h)
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
def ru_bi(vd,vd0,vd1,vd2,omegatot,d1_h):
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
            vdmin, mdmin = vminternal(ec1min,ec2min,omegatot,d1_h)-vd
            ec2mid = -eud-rmid*(1-d1_h)
            ec1mid = -eud+rmid*d1_h
            vdmid, mdmid = vminternal(ec1mid,ec2mid,omegatot,d1_h)-vd
            if abs(vdmid)<error:
                break
            if vdmin*vdmid>0:
                rmin = rmid
            else:
                rmax = rmid
        ru = abs(rmid)
    elif vd<vd2:
        rmin = (-eud-ecu2)/(1-d1_h)
        rmax = 0
        for _ in range(iterations):
            rmid = (rmin+rmax)/2
            ec2min = ecu2
            ec1min = ec2min+rmin
            vdmin, mdmin = vminternal(ec1min,ec2min,omegatot,d1_h)-vd
            ec2mid = ecu2
            ec1mid = ec2mid+rmid
            vdmid, mdmid = vminternal(ec1mid,ec2mid,omegatot,d1_h)-vd
            if abs(vdmid)<error:
                break
            if vdmin*vdmid>0:
                rmin = rmid
            else:
                rmax = rmid
        ru = abs(rmid)
    return ru
d1_h_list = [0.05,0.10,0.15,0.20]
lamda_list = [0,10,20,30,40]
omegatot_list = numpy.arange(0,2.05,0.05)
c_1r = 10.0
bars = 8
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
    data[d1_h] = {}
    for lamda in lamda_list:
        data[d1_h][lamda] = {'curves': {}}
        for omegatot in omegatot_list:
            vd_list = []
            md_max_list = []
            vd0 = -omegatot
            r1 = (-eud-ecu2)/(1-d1_h)
            ec21 = ecu2
            ec11 = ec21+r1
            vd1, md1 = vminternal(ec11,ec21,omegatot,d1_h)
            ec12 = ec22 = 3.5
            vd2, md2 = vminternal(ec12,ec22,omegatot,d1_h)
            vd_array = numpy.linspace(vd0,vd2,1500)
            mphi_data = []
            vd_data = []
            for vd in vd_array:
                ru = ru_bi(vd,vd0,vd1,vd2,omegatot,d1_h)
                m_1st_max = 0
                if ru>1e-6:
                    step = ru/600
                    for Phi in numpy.arange(step,ru+step/2,step):
                        ximid, minternalmid = ximinternalmid(Phi,vd,omegatot,d1_h)
                        if minternalmid is not None:
                            if vd>0:
                                m_2nd = vd*lamda**2/c_1r*Phi/1000
                            else:
                                m_2nd = 0
                            m_1st = minternalmid-m_2nd
                            if m_1st>m_1st_max:
                                m_1st_max=m_1st
                vd_list.append(vd)
                md_max_list.append(m_1st_max)
            data[d1_h][lamda]['curves'][omegatot] = {'vd': vd_list, 'md': md_max_list}
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
for d1_h in d1_h_list:
    for lamda in lamda_list:
        plotdata = data[d1_h][lamda]['curves']
        pyplot.figure(figsize=(11.2,6.6))
        pyplot.xlim(0,1.2)
        pyplot.ylim(0,3)
        ax = pyplot.gca()
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.figure.canvas.draw()
        for omegatot, vd_md in plotdata.items():
            if len(vd_md['md']) == 0: continue
            omega_int = int(round(omegatot*100))
            is_main_curve = (omega_int%40==0) 
            if is_main_curve:
                pyplot.plot(vd_md['md'],vd_md['vd'],linewidth=1,color='k',linestyle='-')
                max_md = max(vd_md['md'])
                max_idx = vd_md['md'].index(max_md)
                vd_nose = vd_md['vd'][max_idx]
                top_idx = max_idx
                while top_idx < len(vd_md['md']) - 1:
                    if vd_md['md'][top_idx+1] < 1e-3:
                        break
                    top_idx += 1
                vd_top = vd_md['vd'][top_idx]
                label_step = omega_int//40
                is_alternate = (label_step%2!=0)
                placement_ratio = 0.60 if is_alternate else 0.40
                target_vd = vd_nose+placement_ratio*(vd_top-vd_nose)
                min_diff = float('inf')
                idx = max_idx
                for i in range(max_idx,top_idx+1):
                    diff = abs(vd_md['vd'][i]-target_vd)
                    if diff < min_diff:
                        min_diff = diff
                        idx = i
                x_text, y_text = vd_md['md'][idx], vd_md['vd'][idx]
                idx_prev = max(max_idx, idx - 15)
                idx_next = min(top_idx, idx + 15)
                if idx_prev >= idx_next:
                    angle = 0
                else:
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
                if omega_int == 200:
                    label_text = '$\omega_{tot}=2.0$'
                else:
                    label_text = '$\omega_{tot}$'+f'={omegatot:.1f}'
                y_offset = 0.01
                pyplot.text(x_text,y_text+y_offset,label_text,
                            fontsize=7,
                            color='k',
                            rotation=angle,
                            rotation_mode='anchor',
                            verticalalignment='bottom',
                            horizontalalignment='center')
            else:
                pyplot.plot(vd_md['md'],vd_md['vd'],linewidth=0.5,color='silver',linestyle='-',alpha=0.8)
        pyplot.grid(True, which='major', linestyle='-', linewidth=0.5, color="#9c9c9c")
        pyplot.savefig(f'plot_d1h_{d1_h}_lamda_{lamda}.png', format='png', dpi=300, bbox_inches='tight')
        pyplot.close()