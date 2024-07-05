-- One temperature reaction rate for Saturn entry condition
-- Author: Yu Liu
-- Date: 26/04/2018
-- Unit: cgs

-- Ref: Direct Simulation Monte Carlo Shock Simulation of Saturn Entry Probe Conditions --Kyle J. Higdon and Brett A. Cruden, 2018
--

-- RJG, 2024-01-15
-- Let's limit the lower end for computing reactions.
-- The reason is that the rate evaluation for the last
-- reaction: He + e- <=> He+ + e- + e-
-- is poorly behaved numerically (causing underflow in exp())
-- at temperatures below 400 K.
-- This is telling us that the forward rate is extremely slow
-- and the backward rate is extremely fast.
-- In other words, it's hard to ionize Helium at low temperatures.
Config{
    tempLimits = {lower=500.0, upper=50000.0}
}

-- H2 dissociation
Reaction{
        'H2 + He <=> H + H + He',
        fr={'Arrhenius', A=4.16963e+18, n=-1.0, C=51958.0465},
        label='r7'
}

Reaction{
        'H2 + H2 <=> H + H + H2',
        fr={'Arrhenius', A=1.04e+19, n=-1.0, C=51958.0465},
        label='r8'
}

Reaction{
        'H2 + H <=> H + H + H',
        fr={'Arrhenius', A=8.34649e+19, n=-1.0, C=51958.0465},
        label='r9'
}

Reaction{
        'H2 + H+ <=> H + H + H+',
        fr={'Arrhenius', A=8.34649e+19, n=-1.0, C=51958.0465},
        label='r10'
}

Reaction{
        'H2 + e- <=> H + H + e-',
        fr={'Arrhenius', A=8.34649e+19, n=-1.0, C=51958.0465},
        label='r11'
}

Reaction{
        'H2 + He+ <=> H + H + He+',
        fr={'Arrhenius', A=4.16963e+18, n=-1.0, C=51958.0465},
        label='r12'
}
--  Ionization

Reaction{
        'H + e- <=> H+ + e- + e-',
        fr={'Arrhenius', A=57.0585e+13, n=0.5, C=157799.7563},
        label='r1'
}

--reaction{
--      'H + e- <=> H+ + e- + e-',
--      fr={'Arrhenius', A=4.11303e+13, n=0.5, T_a=116099.7877},        
--      label='r3'
--}
Reaction{
        'H + H <=> H+ + e- + H',
        fr={'Arrhenius', A=1.541632e+12, n=0.5, C=116099.7877},
        label='r5'
}

Reaction{
        'H + He <=> H+ + e- + He',
        fr={'Arrhenius', A=1.21990665e+12, n=0.5, C=116099.7877},
        label='r6'
}

Reaction{
        'He + e- <=> He+ + e- + e-',
        fr={'Arrhenius', A=33.2715e+13, n=0.5, C=285299.9835},
        label='r2'
}

--reaction{
--      'He + e- <=> He+ + e- + e-',
--      fr={'Arrhenius', A=2.24018e+13, n=0.5, T_a=232100.3466},        
--      label='r4'
--}

