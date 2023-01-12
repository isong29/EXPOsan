# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, numpy as np, qsdsan as qs
from qsdsan import System
from qsdsan.sanunits import ADMtoASM, ASMtoADM
# from exposan.bsm1 import create_system as create_bsm1_system
from exposan.interface import results_path
from exposan.adm import default_init_conds

from exposan.bsm1 import data_path


from qsdsan import processes as pc
Q = 18446           # influent flowrate [m3/d]
Temp = 273.15+20    # temperature [K]

V_an = 1000    # anoxic zone tank volume
V_ae = 1333    # aerated zone tank volume
Q_was = 385    # sludge wastage flowrate
Q_ras = 18446    # recycle sludge flowrate
biomass_IDs = ('X_BH', 'X_BA')

WasteStream = qs.WasteStream

default_inf_kwargs = {
    'concentrations': {
        'S_S':69.5,
        'X_BH':28.17,
        'X_S':202.32,
        'X_I':51.2,
        'S_NH':31.56,
        'S_I':30,
        'S_ND':6.95,
        'X_ND':10.59,
        'S_ALK':7*12,
        },
    'units': ('m3/d', 'mg/L'),
    }

default_asm_kwargs = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3,
    eta_g=0.8, eta_h=0.8, k_h=3.0, K_X=0.1, mu_A=0.5,
    K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05, fr_SS_COD=0.75,
    path=os.path.join(data_path, '_asm1.tsv'),
    )



from qsdsan.utils import load_data
def batch_init(sys, path, sheet):
    df = load_data(path, sheet)
    dct = df.to_dict('index')
    u = sys.flowsheet.unit # unit registry
    for k in [u.A1, u.A2, u.O1, u.O2, u.O3]:
        k.set_init_conc(**dct[k._ID])
    c1s = {k:v for k,v in dct['C1_s'].items() if v>0}
    c1x = {k:v for k,v in dct['C1_x'].items() if v>0}
    tss = [v for v in dct['C1_tss'].values() if v>0]
    u.C1.set_init_solubles(**c1s)
    u.C1.set_init_sludge_solids(**c1x)
    u.C1.set_init_TSS(tss)

su = qs.sanunits

__all__ = ('create_system',)

# fermenters = ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro')
# methanogens = ('X_ac', 'X_h2')
# biomass_IDs = (*fermenters, *methanogens)

def create_system(flowsheet=None):
    flowsheet = flowsheet or qs.Flowsheet('interface')
    qs.main_flowsheet.set_flowsheet(flowsheet)

    # Components and stream
    pc.create_asm1_cmps()
    
    wastewater = WasteStream('wastewater', T=Temp)
    inf_kwargs = default_inf_kwargs
    wastewater.set_flow_by_concentration(Q, **inf_kwargs)
    
    effluent = WasteStream('effluent', T=Temp)
    WAS = WasteStream('WAS', T=Temp)
    RWW = WasteStream('RWW', T=Temp)
    RAS = WasteStream('RAS', T=Temp)
    
    # Process models

    aer1 = aer2 = pc.DiffusedAeration('aer1', 'S_O', KLa=240, DOsat=8.0, V=V_ae)
    aer3 = pc.DiffusedAeration('aer3', 'S_O', KLa=84, DOsat=8.0, V=V_ae)
    asm_kwargs = default_asm_kwargs
    asm1 = pc.ASM1(**asm_kwargs)
    
    # Create unit operations
    A1 = su.CSTR('A1', ins=[wastewater, RWW, RAS, ''], V_max=V_an,
                  aeration=None, suspended_growth_model=asm1)
    
    A2 = su.CSTR('A2', A1-0, V_max=V_an,
                  aeration=None, suspended_growth_model=asm1)
    
    O1 = su.CSTR('O1', A2-0, V_max=V_ae, aeration=aer1,
                  DO_ID='S_O', suspended_growth_model=asm1)
    
    O2 = su.CSTR('O2', O1-0, V_max=V_ae, aeration=aer2,
                  DO_ID='S_O', suspended_growth_model=asm1)
    
    O3 = su.CSTR('O3', O2-0, [RWW, 'treated'], split=[0.6, 0.4],
                 V_max=V_ae, aeration=aer3,
                  DO_ID='S_O', suspended_growth_model=asm1)
    
    C1 = su.FlatBottomCircularClarifier('C1', O3-1, [effluent, RAS, WAS],
                                        underflow=Q_ras, wastage=Q_was, surface_area=1500,
                                        height=4, N_layer=10, feed_layer=5,
                                        X_threshold=3000, v_max=474, v_max_practical=250,
                                        rh=5.76e-4, rp=2.86e-3, fns=2.28e-3)



    # bsm1_sys = create_bsm1_system(flowsheet=flowsheet)
    unit = flowsheet.unit
    stream = flowsheet.stream

    thermo_asm1 = qs.get_thermo() # ASM1 components loaded by the bsm1 module
    cmps_asm1 = thermo_asm1.chemicals
    
    # Subsequent units should be using ADM1 components
    cmps_adm1 = qs.processes.create_adm1_cmps()
    thermo_adm1 = qs.get_thermo()
    adm1 = qs.processes.ADM1()
    cmps_adm1.X_I.i_N = cmps_asm1.X_I.i_N    
    
    J1 = ASMtoADM('J1', upstream=stream.WAS, thermo=thermo_adm1, isdynamic=True, adm1_model=adm1) # WAS is C1.outs[2]
    AD1 = qs.sanunits.AnaerobicCSTR('AD1', ins=J1.outs[0], outs=('biogas', 'ad_eff'), isdynamic=True ,model=adm1,                                    
                                    retain_cmps=[i for i in cmps_adm1.IDs if i.startswith('X_')])
    AD1.set_init_conc(**default_init_conds)
    AD1.T = Temp
    J2 = ADMtoASM('J2', upstream=AD1-1, thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
    
    # Subsequent units should be using ASM1 components
    qs.set_thermo(thermo_asm1)
    # stream.RWW.disconnect_sink() # disconnect from A1 to avoid replacement warning
    
    #!!! Temporary fix (disconnect J2.outs[0] from M1)
    # as the current design will run into errors for long-term simulation
    # S1 = qs.sanunits.HydraulicDelay('S1', ins=[stream.RWW,], isdynamic=True)
    
    # M1 = qs.sanunits.Mixer('M1', ins=[stream.RWW,], isdynamic=True)
    # M1 = qs.sanunits.Mixer('M1', ins=[stream.RWW, J2.outs[0]], isdynamic=True)
    # M1 = qs.sanunits.Mixer('M1', ins=[stream.wastewater, stream.RWW, 
    #                                    J2.outs[0], 
    #                                   stream.RAS],
    #                        isdynamic=True)
    
    A1 = unit.A1

    
    # breakpoint()
    # A1.ins[0] = M1.outs[0]
    A1.ins[-1] = J2.outs[0]
    # unit.A2.ins[0] = A1.outs[0]
    # unit.A1.ins[1] = M1.outs[0]
    # unit.A1.ins[1] = stream.RWW
    
    sys = flowsheet.create_system('interface_sys')
    batch_init(sys, os.path.join(data_path, 'initial_conditions.xlsx'), 'default')
    

    # sys = System('interface_sys', path=(*bsm1_sys.units, J1, AD1, J2, S1))
    # sys = System('interface_sys', path=(*bsm1_sys.units, J1, AD1, J2, S1),
    #              recycle=(S1-0, stream.RAS),
    #              )
    
    # sys = System('interface_sys', path=(*bsm1_sys.units, J1, AD1, J2, M1))
    # sys = System('interface_sys', path=(*bsm1_sys.units, J1, AD1, J2, M1), 
    #              recycle=(stream.RAS, stream.RWW))
    # sys = System('interface_sys', path=(M1, *bsm1_sys.units, J1, AD1, J2,), 
    #               recycle=(M1-0, stream.RAS))
    
    sys.set_dynamic_tracker(unit.A1, unit.C1, J1, AD1, J2)
    # sys.set_dynamic_tracker(unit.A1, unit.C1, J1, AD1, J2, M1)
    
    
    sys.set_tolerance(rmol=1e-10, mol=1e-5)
    sys.maxiter = 5000
    sys.diagram()
    
    return sys


if __name__ == '__main__':
    t = 50
    t_step = 1
    t_eval = np.arange(0, t+t_step, t_step)
    method = 'RK23'
    sys = create_system()
    sys.simulate(
        state_reset_hook='reset_cache',
        t_span=(0, t),
        # t_eval=t_eval,
        method=method,
        )
    # # Just to test a random state
    # states = ('S_su',)
    # AD1 = sys.flowsheet.unit.AD1
    # fig, ax = AD1.scope.plot_time_series(states)
    
    # # Output all states, #!!! seems to have problems
    sys.scope.export(os.path.join(results_path, f'states_{t}_{method}.xlsx'), 
                      t_eval=t_eval)