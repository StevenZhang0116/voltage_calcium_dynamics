COMMENT
gpulse1.mod, patterned after ipulse1.mod
Generates a train of conductance pulses
with a specified reversal potential, 
as if it were a triggerable membrane conductance.
User specifies duration of pulse, interpulse interval (ton and toff),
and number of pulses.
5/27/2008 NTC
ENDCOMMENT

NEURON {
	POINT_PROCESS Gpulse1
	RANGE del, ton, toff, num, e, gmax, g, i
	NONSPECIFIC_CURRENT i
}

UNITS {
	(mV) = (millivolt)
	(nA) = (nanoamp)
	(uS) = (microsiemens)
}

PARAMETER {
	del (ms)
	ton (ms) <0, 1e9>	: duration of ON phase
	toff (ms) <0, 1e9>	: duration of OFF phase
	num			: how many to deliver
	e (mV)			: reversal potential
	gmax (uS) <0, 1e9>	: how big
}

ASSIGNED {
	v (mV)
	g (uS)		: conductance at present time
:	ival (nA)
	i (nA)
	on
	tally		: how many more to deliver
}

INITIAL {
	g = 0
	i = 0
:	ival = 0
	tally = num
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

BREAKPOINT {
: printf("%g\n", t)
	i = g*(v-e)
}

NET_RECEIVE (w) {
	: ignore any but self-events with flag == 1
	if (flag == 1) {
		if (on == 0) {
			: turn it on
:			ival = amp
			g = gmax
			on = 1
			: prepare to turn it off
			net_send(ton, 1)
		} else {
			: turn it off
:			ival = 0
			g = 0
			on = 0
			if (tally > 0) {
				: prepare to turn it on again
				net_send(toff, 1)
				tally = tally - 1
			}
		}
	}
}
