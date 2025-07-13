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
    metal = "Fe",
    facet = "bcc110",
)
        
entry(
    index = 0,
    label = "C",
    molecule = """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0.0774496,0.0207201,-3.33637e-05,3.47705e-08,-1.39054e-11,-7509.59,19.4592], Tmin=(298.15,'K'), Tmax=(739.823,'K')), NASAPolynomial(coeffs=[2.07402,0.00682633,1.08894e-06,-1.93703e-09,4.11971e-13,-7720.21,11.0017], Tmin=(739.823,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(-61.8729,'kJ/mol'), Cp0=(0.000344693,'J/(mol*K)'), CpInf=(0.00112026,'J/(mol*K)')),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: -13.510609503569043 [kcal/mol]
Sf298: 49.427609990392114 [cal/(mol-K)]
C 10.71809929 10.6520347 10.6520347
H 11.24749643 11.39927707 11.39927707
H 11.43380999 10.06585349 10.06585349
H 10.16128142 9.99634732 9.99634732
H 10.02830394 11.1471665 11.1471665
""",
)


entry(
    index = 1,
    label = "N",
    molecule = """
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.172197,0.0263519,-5.64738e-05,6.07965e-08,-2.4363e-11,-4176.72,19.3927], Tmin=(298.15,'K'), Tmax=(638.576,'K')), NASAPolynomial(coeffs=[3.9545,0.000502195,4.24729e-06,-2.59648e-09,4.55405e-13,-4703.76,1.33472], Tmin=(638.576,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(-34.2551,'kJ/mol'), Cp0=(0.000344693,'J/(mol*K)'), CpInf=(0.000861783,'J/(mol*K)')),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: -6.851321108962932 [kcal/mol]
Sf298: 48.180368724519056 [cal/(mol-K)]
N 10.59664039 10.77282591 10.77282591
H 11.52795165 10.66565161 10.66565161
H 10.03350828 10.01748155 10.01748155
H 10.21945027 11.63595359 11.63595359
""",
)


entry(
    index = 2,
    label = "N#[Pt]",
    molecule = 
"""
1 X u0 p0 c0 {2,T}
2 N u0 p1 c0 {1,T}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0.628962,0.0111444,-2.178e-05,1.9796e-08,-6.89307e-12,1.71584e+07,-3.93665], Tmin=(298.15,'K'), Tmax=(703.852,'K')), NASAPolynomial(coeffs=[2.29041,0.00170239,-1.65772e-06,7.36799e-10,-1.23445e-13,1.71582e+07,-11.3687], Tmin=(703.852,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(142664,'kJ/mol'), Cp0=(2.86452e-11,'J/(mol*K)'), CpInf=(24.9434,'J/(mol*K)'), comment="""metal = Fe,facet = bcc110"""),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 34098.317133757104 [kcal/mol]
Sf298: 4.295567198351038 [cal/(mol-K)]
28
Fe 0.0 0.0 8.0
Fe 2.831 0.0 8.0
Fe 5.662 0.0 8.0
Fe 1.4155 2.0018193 8.0
Fe 4.2465 2.0018193 8.0
Fe 7.0775 2.0018193 8.0
Fe 2.831 4.0036386 8.0
Fe 5.662 4.0036386 8.0
Fe 8.493 4.0036386 8.0
Fe 0.0 2.0018193 10.0018193
Fe 2.831 2.0018193 10.0018193
Fe 5.662 2.0018193 10.0018193
Fe 1.4155 4.0036386 10.0018193
Fe 4.2465 4.0036386 10.0018193
Fe 7.0775 4.0036386 10.0018193
Fe 2.831 6.00545789 10.0018193
Fe 5.662 6.00545789 10.0018193
Fe 8.493 6.00545789 10.0018193
Fe -0.13759877 0.66177529 11.86142089
Fe 2.83554003 0.45299144 11.84149257
Fe 5.78927281 0.65074188 11.85147313
Fe 1.46281162 2.5369144 11.81973618
Fe 4.19139026 2.53118178 11.82205039
Fe 7.06206306 2.72869964 11.93060005
Fe 2.83694055 4.58864049 11.85862949
Fe 5.33591969 4.66643495 11.94032508
Fe 8.812759 4.65987273 11.9358411
N 7.07660812 4.53126275 12.49138882
""",
)


entry(
    index = 3,
    label = "N=[Pt]",
    molecule = 
"""
1 X u0 p0 c0 {2,D}
2 N u0 p1 c0 {1,D} {3,S}
3 H u0 p0 c0 {2,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.602222,0.0234081,-4.20726e-05,3.60478e-08,-1.17931e-11,1.72997e+07,0.618092], Tmin=(298.15,'K'), Tmax=(755.192,'K')), NASAPolynomial(coeffs=[3.18628,0.00334046,-2.21091e-06,8.5667e-10,-1.4268e-13,1.72991e+07,-16.5953], Tmin=(755.192,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(143838,'kJ/mol'), Cp0=(5.39236e-11,'J/(mol*K)'), CpInf=(49.8868,'J/(mol*K)'), comment="""metal = Fe,facet = bcc110"""),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 34379.09184254115 [kcal/mol]
Sf298: 5.145688820699216 [cal/(mol-K)]
29
Fe 0.0 0.0 8.0
Fe 2.831 0.0 8.0
Fe 5.662 0.0 8.0
Fe 1.4155 2.0018193 8.0
Fe 4.2465 2.0018193 8.0
Fe 7.0775 2.0018193 8.0
Fe 2.831 4.0036386 8.0
Fe 5.662 4.0036386 8.0
Fe 8.493 4.0036386 8.0
Fe 0.0 2.0018193 10.0018193
Fe 2.831 2.0018193 10.0018193
Fe 5.662 2.0018193 10.0018193
Fe 1.4155 4.0036386 10.0018193
Fe 4.2465 4.0036386 10.0018193
Fe 7.0775 4.0036386 10.0018193
Fe 2.831 6.00545789 10.0018193
Fe 5.662 6.00545789 10.0018193
Fe 8.493 6.00545789 10.0018193
Fe 0.0529069 0.52694155 11.82678456
Fe 2.7923538 0.51329705 11.82584007
Fe 4.41090264 -1.33100437 11.87130206
Fe 1.41264464 2.57357503 11.84610593
Fe 3.93208151 2.68174822 11.94033034
Fe 5.65622153 0.74946018 11.92462848
Fe 2.67202044 4.67320101 11.87027352
Fe 5.65880527 4.38687761 11.8334976
Fe 7.38372788 2.68139794 11.93281524
N 5.65796195 2.5314578 12.82381751
H 5.66628862 2.50490022 13.8592355
""",
)


entry(
    index = 4,
    label = "N[Pt]",
    molecule = 
"""
1 X u0 p0 c0 {2,S}
2 N u0 p1 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.789977,0.0280871,-4.76641e-05,4.05643e-08,-1.3328e-11,1.74493e+07,1.48217], Tmin=(298.15,'K'), Tmax=(750.155,'K')), NASAPolynomial(coeffs=[3.38203,0.00584177,-3.18426e-06,1.0363e-09,-1.55143e-13,1.74487e+07,-17.4462], Tmin=(750.155,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(145082,'kJ/mol'), Cp0=(1.07255e-05,'J/(mol*K)'), CpInf=(74.8302,'J/(mol*K)'), comment="""metal = Fe,facet = bcc110"""),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 34676.70987117844 [kcal/mol]
Sf298: 7.087806509776206 [cal/(mol-K)]
30
Fe 0.0 0.0 8.0
Fe 2.831 0.0 8.0
Fe 5.662 0.0 8.0
Fe 1.4155 2.0018193 8.0
Fe 4.2465 2.0018193 8.0
Fe 7.0775 2.0018193 8.0
Fe 2.831 4.0036386 8.0
Fe 5.662 4.0036386 8.0
Fe 8.493 4.0036386 8.0
Fe 0.0 2.0018193 10.0018193
Fe 2.831 2.0018193 10.0018193
Fe 5.662 2.0018193 10.0018193
Fe 1.4155 4.0036386 10.0018193
Fe 4.2465 4.0036386 10.0018193
Fe 7.0775 4.0036386 10.0018193
Fe 2.831 6.00545789 10.0018193
Fe 5.662 6.00545789 10.0018193
Fe 8.493 6.00545789 10.0018193
Fe 0.03741592 -0.48914617 11.79287658
Fe 2.92453378 -0.6139656 11.85604292
Fe 5.54632602 -0.61842095 11.85327239
Fe 1.66696818 1.38343146 11.84399285
Fe 4.2115747 1.43718857 11.85463923
Fe 6.79451448 1.34196791 11.96079996
Fe 2.82388211 3.48373599 11.82196881
Fe 5.58356007 3.48685678 11.8052431
Fe 8.44862271 3.306475 11.95141457
N 7.70412647 2.28109225 13.45161017
H 7.14790741 2.83816613 14.11337081
H 8.35895494 1.7060112 13.99870612
""",
)


entry(
    index = 5,
    label = "O",
    molecule = """
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0.246628,0.0254069,-6.19234e-05,6.79605e-08,-2.72995e-11,-28351.3,16.6214], Tmin=(298.15,'K'), Tmax=(637.488,'K')), NASAPolynomial(coeffs=[4.83253,-0.00336674,5.77771e-06,-2.8364e-09,4.63376e-13,-28936,-3.43848], Tmin=(637.488,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(-234.937,'kJ/mol'), Cp0=(0.000344693,'J/(mol*K)'), CpInf=(0.000603199,'J/(mol*K)')),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: -54.79674419465377 [kcal/mol]
Sf298: 46.48839325422904 [cal/(mol-K)]
O 10.81193598 10.57360488 10.57360488
H 10.0496051 9.98800591 9.98800591
H 11.57540752 9.98883554 9.98883554
""",
)


entry(
    index = 6,
    label = "O=O",
    molecule = """
1 O u0 p2 c0 {2,D}
2 O u0 p2 c0 {1,D}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0.569878,0.0204858,-5.25591e-05,6.19645e-08,-2.72482e-11,15675.7,16.7299], Tmin=(298.15,'K'), Tmax=(575.298,'K')), NASAPolynomial(coeffs=[3.57979,-0.000444418,2.01999e-06,-1.29061e-09,2.4312e-13,15329.5,3.87321], Tmin=(575.298,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(131.157,'kJ/mol'), Cp0=(0.000301607,'J/(mol*K)'), CpInf=(0.000387797,'J/(mol*K)')),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 32.59198340663045 [kcal/mol]
Sf298: 48.17077782292974 [cal/(mol-K)]
O 11.16872454 10.0 10.0
O 9.97220501 10.0 10.0
""",
)


entry(
    index = 7,
    label = "O=[Pt]",
    molecule = 
"""
1 X u0 p0 c0 {2,D}
2 O u0 p2 c0 {1,D}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0.623905,0.0112489,-2.20916e-05,2.01527e-08,-7.03845e-12,2.36168e+07,-3.98389], Tmin=(298.15,'K'), Tmax=(702.191,'K')), NASAPolynomial(coeffs=[2.30553,0.00166966,-1.62856e-06,7.24908e-10,-1.2162e-13,2.36165e+07,-11.5022], Tmin=(702.191,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(196362,'kJ/mol'), Cp0=(6.74236e-12,'J/(mol*K)'), CpInf=(24.9434,'J/(mol*K)'), comment="""metal = Fe,facet = bcc110"""),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 46932.461824273734 [kcal/mol]
Sf298: 4.184512048258376 [cal/(mol-K)]
28
Fe 0.0 0.0 8.0
Fe 2.831 0.0 8.0
Fe 5.662 0.0 8.0
Fe 1.4155 2.0018193 8.0
Fe 4.2465 2.0018193 8.0
Fe 7.0775 2.0018193 8.0
Fe 2.831 4.0036386 8.0
Fe 5.662 4.0036386 8.0
Fe 8.493 4.0036386 8.0
Fe 0.0 2.0018193 10.0018193
Fe 2.831 2.0018193 10.0018193
Fe 5.662 2.0018193 10.0018193
Fe 1.4155 4.0036386 10.0018193
Fe 4.2465 4.0036386 10.0018193
Fe 7.0775 4.0036386 10.0018193
Fe 2.831 6.00545789 10.0018193
Fe 5.662 6.00545789 10.0018193
Fe 8.493 6.00545789 10.0018193
Fe 0.16349241 0.6102666 11.84306957
Fe 2.83400258 0.57092503 11.81790457
Fe 5.46694845 0.65556782 11.86044651
Fe 1.56718638 2.59366754 11.85310585
Fe 4.05467022 2.61797855 11.85723735
Fe 7.09693754 2.43109276 11.86461279
Fe 2.77975587 4.70353624 11.88412068
Fe 5.54574431 4.57375622 11.92993881
Fe 8.6378554 4.5652272 11.9080707
O 7.09334074 4.01354239 12.80253345
""",
)


entry(
    index = 8,
    label = "O[Pt]",
    molecule = 
"""
1 X u0 p0 c0 {2,S}
2 O u0 p2 c0 {1,S} {3,S}
3 H u0 p0 c0 {2,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.0159245,0.0232638,-4.56991e-05,4.20341e-08,-1.46509e-11,2.37658e+07,-1.73076], Tmin=(298.15,'K'), Tmax=(712.509,'K')), NASAPolynomial(coeffs=[3.7342,0.00221058,-1.37697e-06,5.63426e-10,-9.98682e-14,2.37652e+07,-18.5518], Tmin=(712.509,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(197600,'kJ/mol'), Cp0=(8.23524e-09,'J/(mol*K)'), CpInf=(49.8868,'J/(mol*K)'), comment="""metal = Fe,facet = bcc110"""),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 47228.87438215474 [kcal/mol]
Sf298: 6.803983740839437 [cal/(mol-K)]
29
Fe 0.0 0.0 8.0
Fe 2.831 0.0 8.0
Fe 5.662 0.0 8.0
Fe 1.4155 2.0018193 8.0
Fe 4.2465 2.0018193 8.0
Fe 7.0775 2.0018193 8.0
Fe 2.831 4.0036386 8.0
Fe 5.662 4.0036386 8.0
Fe 8.493 4.0036386 8.0
Fe 0.0 2.0018193 10.0018193
Fe 2.831 2.0018193 10.0018193
Fe 5.662 2.0018193 10.0018193
Fe 1.4155 4.0036386 10.0018193
Fe 4.2465 4.0036386 10.0018193
Fe 7.0775 4.0036386 10.0018193
Fe 2.831 6.00545789 10.0018193
Fe 5.662 6.00545789 10.0018193
Fe 8.493 6.00545789 10.0018193
Fe 0.12726763 0.56517577 11.84049183
Fe 2.84603649 0.5787006 11.8130772
Fe 5.44068274 0.67079957 11.88013672
Fe 1.51521491 2.57650792 11.83905801
Fe 4.05513961 2.63868269 11.86595073
Fe 7.08804888 2.40112033 11.85954878
Fe 2.6866077 4.69551046 11.8898211
Fe 5.56870126 4.57582606 11.89265766
Fe 8.63823459 4.57932554 11.9027425
O 7.09325603 3.91070387 13.10015111
H 7.09607635 3.75271843 14.06388339
""",
)


entry(
    index = 9,
    label = "[H][H]",
    molecule = """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[0.408382,0.0214959,-5.33395e-05,5.67745e-08,-2.2102e-11,256.761,9.54822], Tmin=(298.15,'K'), Tmax=(654.873,'K')), NASAPolynomial(coeffs=[4.52538,-0.00365215,4.26524e-06,-1.87051e-09,2.86925e-13,-282.435,-8.57092], Tmin=(654.873,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(2.97588,'kJ/mol'), Cp0=(0.000301607,'J/(mol*K)'), CpInf=(0.000387775,'J/(mol*K)')),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 1.9155568841383364 [kcal/mol]
Sf298: 32.529449943673285 [cal/(mol-K)]
H 10.71973162 10.0 10.0
H 9.97826883 10.0 10.0
""",
)


entry(
    index = 10,
    label = "[Pt]",
    molecule = 
"""
1 X u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
""",
thermo=NASA(polynomials=[NASAPolynomial(coeffs=[-0.807111,0.00885124,-6.97011e-06,9.6513e-10,6.81688e-13,111056,2.72021], Tmin=(298.15,'K'), Tmax=(909.062,'K')), NASAPolynomial(coeffs=[0.138196,0.00662923,-6.5006e-06,2.9653e-09,-5.13132e-13,110804,-2.19054], Tmin=(909.062,'K'), Tmax=(2000,'K'))], Tmin=(298.15,'K'), Tmax=(2000,'K'), E0=(922.851,'kJ/mol'), Cp0=(3.04493e-33,'J/(mol*K)'), CpInf=(24.9434,'J/(mol*K)'), comment="""metal = Fe,facet = bcc110"""),    
longDesc = u"""Calculated using Pynta https://doi.org/10.1021/acs.jcim.3c00948. 
                   Thermochemistry computed using approach from Blondal et. al in https://doi.org/10.1021/acscatal.2c03378. 
If you use this library in your work, please cite the publications mentioned above.
Hf298: 220.87556233491327 [kcal/mol]
Sf298: 0.9142465196911341 [cal/(mol-K)]
28
Fe 0.0 0.0 8.0
Fe 2.831 0.0 8.0
Fe 5.662 0.0 8.0
Fe 1.4155 2.0018193 8.0
Fe 4.2465 2.0018193 8.0
Fe 7.0775 2.0018193 8.0
Fe 2.831 4.0036386 8.0
Fe 5.662 4.0036386 8.0
Fe 8.493 4.0036386 8.0
Fe 0.0 2.0018193 10.0018193
Fe 2.831 2.0018193 10.0018193
Fe 5.662 2.0018193 10.0018193
Fe 1.4155 4.0036386 10.0018193
Fe 4.2465 4.0036386 10.0018193
Fe 7.0775 4.0036386 10.0018193
Fe 2.831 6.00545789 10.0018193
Fe 5.662 6.00545789 10.0018193
Fe 8.493 6.00545789 10.0018193
Fe 0.1045071 0.52291672 11.82009719
Fe 2.87432711 0.60880725 11.84593795
Fe 5.49829054 0.61001753 11.85418538
Fe 1.49506109 2.58468476 11.84349061
Fe 4.15291266 2.65609254 11.86437835
Fe 7.08144435 2.43998261 11.84941559
Fe 2.69932593 4.66854616 11.86631341
Fe 5.8017916 4.55699487 11.86514829
Fe 8.44961791 4.54224879 11.86053038
H 7.11549289 3.82640367 12.87701615
""",
)

