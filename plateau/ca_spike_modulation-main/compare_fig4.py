from Fig4_3_two_variables import *
from Fig4_2_four_variables import *

plt.figure(figsize=(9,9))
plt.plot(t_fig42,v_fig42,label='4 params')
plt.plot(t_fig43,v_fig43,label='2 params')
plt.xlim([-5, 60])
plt.xlabel('Time (ms)',fontsize=16)
plt.ylabel('Vm (mV)',fontsize=16)
plt.legend(fontsize=14)
plt.show()

plt.figure(figsize=(9,9))
plt.plot(t_fig42,Ca_HVA_m_t_fig42,label='m_{Ca_{HVA}} - 4 params')
plt.plot(t_fig42,Im_m_t_fig42,label='n_{I_m} - 4 params')
plt.plot(t_fig43,Ca_HVA_m_t_fig43,'--',label='m_{Ca_{HVA}} - 2 params')
plt.plot(t_fig43,Im_m_t_fig43,'--',label='n_{I_m} - 2 params')
plt.xlim([-5, 60])
plt.xlabel('Time (ms)',fontsize=16)
plt.ylabel('Gate',fontsize=16)
plt.legend(fontsize=14)
plt.show()
