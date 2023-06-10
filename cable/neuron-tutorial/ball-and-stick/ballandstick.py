from neuron import h,gui
from neuron.units import mV, ms
h.load_file('stdrun.hoc')

# Define the general class
class Cell:
    def __init__(self, gid, x, y, z, theta):
        self._gid = gid
        self._setup_morphology()
        self.all = self.soma.wholetree()
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        h.define_shape()
        self._rotate_z(theta)
        self._set_position(x, y, z)
        
        # everything below here in this method is NEW
        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        
        self._ncs = []
        self.soma_v = h.Vector().record(self.soma(0.5)._ref_v)

        
    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)
    
    def _set_position(self, x, y, z):
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(i,
                               x - self.x + sec.x3d(i),
                               y - self.y + sec.y3d(i),
                               z - self.z + sec.z3d(i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
        
    def _rotate_z(self, theta):
        """Rotate the cell about the Z axis."""
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))


class BallAndStick(Cell):
    name = 'BallAndStick'
    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.all = [self.soma,self.dend]
        self.dend.connect(self.soma)
        self.soma.L = self.soma.diam = 12.6157
        self.dend.L = 200
        self.dend.diam = 1

    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100 # Axial resistance
            sec.cm = 1 # Membrane capacitance
        self.soma.insert('hh')
        for seg in self.soma:
            seg.hh.gnabar = 0.12 # Sodium conductance
            seg.hh.gkbar = 0.036 # Potassium conductance
            seg.hh.gl = 0.0003 # Leak Conductance
            seg.hh.el = -54.3 # Reversal potential
        self.dend.insert('pas')
        for seg in self.dend:
            seg.pas.g = 0.001 # Passive conductance
            seg.pas.e = -65 # Leak reversal potential

        self.syn = h.ExpSyn(self.dend(0.5))
        self.syn.tau = 2 * ms