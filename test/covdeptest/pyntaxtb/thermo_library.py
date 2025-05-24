#!/usr/bin/env python
# encoding: utf-8

name = "thermo_library"
shortDesc = u""
longDesc = u""""""
#

entry(
    index = 1,
    label = "vacant",
    molecule =
"""
1 X  u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
            NASAPolynomial(coeffs=[
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00,   0.000000000E+00,
             0.000000000E+00,   0.000000000E+00,   0.000000000E+00], Tmin=(1000.0,'K'), Tmax=(3000.0, 'K')),
        ],
        Tmin = (298.0, 'K'),
        Tmax = (3000.0, 'K'),
    ),
    metal = "Cu",
    facet = "111",
)
        
entry(
    index = 0,
    label = "[H][H]",
    molecule = """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0.380556,0.0214452,-5.25131e-05,5.51192e-08,-2.1187e-11,259.225,9.76209], Tmin=(298.15,'K'), Tmax=(662.397,'K')), NASAPolynomial(coeffs=[4.50766,-0.0034782,3.92854e-06,-1.68886e-09,2.54334e-13,-287.506,-8.44869], Tmin=(662.397,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(2.97996,'kJ/mol'), Cp0=(0.000301607,'J/(mol*K)'), CpInf=(0.000387786,'J/(mol*K)')),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 1.908357965270702 [kcal/mol]
Sf298: 32.65690066676592 [cal/(mol-K)]
H 10.7318129 10.0 10.0
H 9.96618755 10.0 10.0
""",
)

