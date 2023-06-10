import time
import pickle
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import phys_funs as pf
import channel_funs as cf
import dynamics_funs as df
from cfg import *
from constants import *


# =========== parameters ===========
file = open("nexus_point_neuron_parameters.pkl",'rb')
nexus_params = pickle.load(file)
file.close()

general_params = cfg_fig_4_2["h"]
tmin = general_params["tmin"]
tmax = general_params["tmax"]
dt = general_params["dt"]
T = np.arange(tmin, tmax, dt)

phys_params = cfg_fig_4_2["phys"]
E_K = phys_params["E_K"]
E_Ca = phys_params["E_Ca"]
E_l = phys_params["E_l"]
Cm = phys_params["Cm"]
A = phys_params["A"]

sim_params = cfg_fig_4_2["sim"]
amps = sim_params["amps"]
mul_Cas = sim_params["mul_Cas"]
mul_Ks = sim_params["mul_Ks"]
ts = sim_params["ts"]
tau1 = sim_params["tau1"]
tau2 = sim_params["tau2"]
amp_syn = sim_params["amp_syn"]
stim_delay = sim_params["stim_delay"]

# ============ run simulation =============
dic_log = {}
for amp in amps:
    for mul_Ca in mul_Cas:
        for mul_K in mul_Ks:
            tic = time.time()
            data = {}

            gCa_HVAbar = nexus_params['gCa_HVAbar_Ca_HVA'] * mul_Ca * TO_MS
            gImbar = nexus_params['gImbar_Im'] * mul_K * TO_MS
            gl_pas = nexus_params['g_pas'] * TO_MS

            # State vector parameters: v, all gates
            V_init = np.copy(E_l)
            Y = np.array([V_init, cf.Ca_HVA_m_inf(V_init), cf.Ca_HVA_h_inf(V_init), cf.Im_m_inf(V_init)])

            # Solve ODE system
            Vy = odeint(df.compute_derivatives_four_variables_model, Y, T, args=(amp, amp_syn, ts, tau1, tau2,
                                                                                 gCa_HVAbar, gImbar, gl_pas,
                                                                                 E_Ca, E_K, E_l, A, Cm))

            # get derivatives
            dy_list = []
            for i in range(len(T)):
                dy = df.compute_derivatives_four_variables_model(Vy[i], T[i], amp, amp_syn, ts, tau1, tau2,
                                                                 gCa_HVAbar, gImbar, gl_pas,
                                                                 E_Ca, E_K, E_l, A, Cm)
                dy_list.append(dy)
            Vdy = np.array(dy_list)
            v_dot_t = Vdy[:, 0]
            Ca_HVA_m_dot_t = Vdy[:, 1]
            Ca_HVA_h_dot_t = Vdy[:, 2]
            Im_m_dot_t = Vdy[:, 3]

            v_t = Vy[:, 0]
            Ca_HVA_m_t = Vy[:, 1]
            Ca_HVA_h_t = Vy[:, 2]
            Im_m_t = Vy[:, 3]

            gCa_HVA_t = pf.calc_conductance('HVA', gCa_HVAbar, [Ca_HVA_m_t, Ca_HVA_h_t])
            gIm_t = pf.calc_conductance('Im', gImbar, [Im_m_t])
            gl_pas_t = gl_pas * np.ones(np.size(v_t))

            I_Ca_HVA_t = pf.calc_current(v_t, gCa_HVA_t, E_Ca)
            I_Im_t = pf.calc_current(v_t, gIm_t, E_K)
            I_l_t = pf.calc_current(v_t, gl_pas_t, E_l)

            # Input stimulus
            I_ext_per_area = pf.I_ext_constant(amp) / A
            I_ext_syn_vec = np.array([pf.I_ext_syn(t_curr, amp_syn, ts, tau1, tau2) / A for t_curr in T])
            I_ext_sum_vec = I_ext_per_area * np.ones_like(I_ext_syn_vec) + I_ext_syn_vec

            data['v'] = v_t
            data['t'] = T

            data['I_ext_per_area'] = I_ext_per_area
            data['I_ext_syn_vec'] = I_ext_syn_vec
            data['I_ext_sum_vec'] = I_ext_sum_vec
            data['I_l_t'] = I_l_t
            data['I_Ca_HVA_t'] = I_Ca_HVA_t
            data['I_Im_t'] = I_Im_t

            data['Ca_HVA_m_t'] = Ca_HVA_m_t
            data['Ca_HVA_h_t'] = Ca_HVA_h_t
            data['Im_m_t'] = Im_m_t

            data['gCa_HVA_t'] = gCa_HVA_t
            data['gIm_t'] = gIm_t
            data['gl_pas_t'] = gl_pas_t

            dic_log[(amp, E_l, mul_Ca, mul_K)] = data

            toc = time.time()
            print('finished %f amp %f E_l %f mul_Ca %f mul_K in %f secs' % (amp, E_l, mul_Ca, mul_K, toc - tic))

# ======== plot ===========
from matplotlib import cm

evenly_spaced_interval = np.linspace(0, 1, 7)
colors = [cm.coolwarm(k) for k in evenly_spaced_interval]

fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(10,10), sharex=True, gridspec_kw={'height_ratios': [1, 1, 0.2]})
for amp in amps:
    for mul_Ca in mul_Cas:
        for mul_K in mul_Ks:
            t = dic_log[(amp, E_l, mul_Ca, mul_K)]['t']
            t_shifted = t - stim_delay

            v = dic_log[(amp, E_l, mul_Ca, mul_K)]['v']
            Ca_HVA_m_t = dic_log[(amp, E_l, mul_Ca, mul_K)]['Ca_HVA_m_t']
            Ca_HVA_h_t = dic_log[(amp, E_l, mul_Ca, mul_K)]['Ca_HVA_h_t']
            Im_m_t1 = dic_log[(amp, E_l, mul_Ca, mul_K)]['Im_m_t']
            I_ext_sum_vec = dic_log[(amp, E_l, mul_Ca, mul_K)]['I_ext_sum_vec']

            ax0.plot(t_shifted, v, color='r')
            line_Ca_HVA, = ax1.plot(t_shifted, Ca_HVA_m_t, label='mCa', color=colors[6])
            ax1.plot(t_shifted, Ca_HVA_h_t, '--', label='hCa', color=colors[6], linestyle='dashed')
            line_Im, = ax1.plot(t_shifted, Im_m_t1, label='mIm', color=colors[0])
            ax2.plot(t_shifted, I_ext_sum_vec, color='r')

ax0.set_ylabel('Vm (mV)', fontsize=16)
ax0.set_xlim([-5, 60])
ax0.set_ylim([-90, 50])

ax1.set_ylabel('gate', fontsize=16)
ax1.set_xlim([-5, 60])
ax1.legend(loc='upper right')

ax2.set_ylabel('I', fontsize=16)
ax2.set_xlabel('ms', fontsize=16)
ax2.set_xlim([-5, 60])

plt.show()
if general_params["save_figs"]:
    fig.savefig("Fig4/fig4_CaHVA_Im_two_channels_model.svg", format='svg')

v_fig42 = v
t_fig42 = t_shifted
Ca_HVA_m_t_fig42 = Ca_HVA_m_t
Im_m_t_fig42 = Im_m_t1
