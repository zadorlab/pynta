import unittest
import os
from pynta.mol import *

class UtilsTest(unittest.TestCase):

    def test_get_full_mol(self):
        test_out_mol_adjlist = """1    X u0 p0 c0 s"ontop" m"terrace" {2,S} {5,S} {8,S} {10,S} {11,S} {13,S}
2    X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {4,S} {14,S}
3    X u0 p0 c0 s"fcc" m"terrace" {2,S} {8,S} {18,S}
4    X u0 p0 c0 s"hcp" m"terrace" {2,S} {13,S} {21,S}
5    X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {7,S} {24,S}
6    X u0 p0 c0 s"fcc" m"terrace" {5,S} {10,S} {27,S}
7    X u0 p0 c0 s"hcp" m"terrace" {5,S} {11,S} {28,S}
8    X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {9,S} {31,S}
9    X u0 p0 c0 s"hcp" m"terrace" {8,S} {10,S} {34,S}
10   X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {9,S} {45,S}
11   X u0 p0 c0 s"bridge" m"terrace" {1,S} {7,S} {12,S} {49,S}
12   X u0 p0 c0 s"fcc" m"terrace" {11,S} {13,S} {50,S}
13   X u0 p0 c0 s"bridge" m"terrace" {1,S} {4,S} {12,S} {52,S}
14   X u0 p0 c0 s"ontop" m"terrace" {2,S} {15,S} {18,S} {20,S} {21,S} {23,S}
15   X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {17,S} {24,S}
16   X u0 p0 c0 s"fcc" m"terrace" {15,S} {20,S} {25,S}
17   X u0 p0 c0 s"hcp" m"terrace" {15,S} {23,S} {30,S}
18   X u0 p0 c0 s"bridge" m"terrace" {3,S} {14,S} {19,S} {31,S}
19   X u0 p0 c0 s"hcp" m"terrace" {18,S} {20,S} {32,S}
20   X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {19,S} {39,S}
21   X u0 p0 c0 s"bridge" m"terrace" {4,S} {14,S} {22,S} {52,S}
22   X u0 p0 c0 s"fcc" m"terrace" {21,S} {23,S} {53,S}
23   X u0 p0 c0 s"bridge" m"terrace" {14,S} {17,S} {22,S} {54,S}
24   X u0 p0 c0 s"ontop" m"terrace" {5,S} {15,S} {25,S} {27,S} {28,S} {30,S}
25   X u0 p0 c0 s"bridge" m"terrace" {16,S} {24,S} {26,S} {39,S}
26   X u0 p0 c0 s"hcp" m"terrace" {25,S} {27,S} {40,S}
27   X u0 p0 c0 s"bridge" m"terrace" {6,S} {24,S} {26,S} {45,S}
28   X u0 p0 c0 s"bridge" m"terrace" {7,S} {24,S} {29,S} {49,S}
29   X u0 p0 c0 s"fcc" m"terrace" {28,S} {30,S} {51,S}
30   X u0 p0 c0 s"bridge" m"terrace" {17,S} {24,S} {29,S} {54,S}
31   X u0 p0 c0 s"ontop" m"terrace" {8,S} {18,S} {32,S} {34,S} {36,S} {38,S}
32   X u0 p0 c0 s"bridge" m"terrace" {19,S} {31,S} {33,S} {39,S}
33   X u0 p0 c0 s"fcc" m"terrace" {32,S} {36,S} {42,S} {55,D}
34   X u0 p0 c0 s"bridge" m"terrace" {9,S} {31,S} {35,S} {45,S}
35   X u0 p0 c0 s"fcc" m"terrace" {34,S} {38,S} {48,S}
36   X u0 p0 c0 s"bridge" m"terrace" {31,S} {33,S} {37,S} {49,S}
37 * X u0 p0 c0 s"hcp" m"terrace" {36,S} {38,S} {51,S} {56,R}
38   X u0 p0 c0 s"bridge" m"terrace" {31,S} {35,S} {37,S} {54,S}
39   X u0 p0 c0 s"ontop" m"terrace" {20,S} {25,S} {32,S} {40,S} {42,S} {44,S}
40   X u0 p0 c0 s"bridge" m"terrace" {26,S} {39,S} {41,S} {45,S}
41   X u0 p0 c0 s"fcc" m"terrace" {40,S} {44,S} {46,S}
42   X u0 p0 c0 s"bridge" m"terrace" {33,S} {39,S} {43,S} {49,S}
43   X u0 p0 c0 s"hcp" m"terrace" {42,S} {44,S} {50,S}
44   X u0 p0 c0 s"bridge" m"terrace" {39,S} {41,S} {43,S} {52,S}
45   X u0 p0 c0 s"ontop" m"terrace" {10,S} {27,S} {34,S} {40,S} {46,S} {48,S}
46   X u0 p0 c0 s"bridge" m"terrace" {41,S} {45,S} {47,S} {52,S}
47   X u0 p0 c0 s"hcp" m"terrace" {46,S} {48,S} {53,S}
48   X u0 p0 c0 s"bridge" m"terrace" {35,S} {45,S} {47,S} {54,S}
49   X u0 p0 c0 s"ontop" m"terrace" {11,S} {28,S} {36,S} {42,S} {50,S} {51,S}
50   X u0 p0 c0 s"bridge" m"terrace" {12,S} {43,S} {49,S} {52,S}
51   X u0 p0 c0 s"bridge" m"terrace" {29,S} {37,S} {49,S} {54,S}
52   X u0 p0 c0 s"ontop" m"terrace" {13,S} {21,S} {44,S} {46,S} {50,S} {53,S}
53   X u0 p0 c0 s"bridge" m"terrace" {22,S} {47,S} {52,S} {54,S}
54   X u0 p0 c0 s"ontop" m"terrace" {23,S} {30,S} {38,S} {48,S} {51,S} {53,S}
55   O u0 p2 c0 {33,D} {56,R}
56 * H u0 p0 c0 {37,R} {55,R}"""
        test_out = Molecule().from_adjacency_list(test_out_mol_adjlist,check_consistency=False)
        slab_struct = """1    X u0 p0 c0 s"ontop" m"terrace" {2,S} {5,S} {8,S} {10,S} {11,S} {13,S}
2    X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {4,S} {14,S}
3    X u0 p0 c0 s"fcc" m"terrace" {2,S} {8,S} {18,S}
4    X u0 p0 c0 s"hcp" m"terrace" {2,S} {13,S} {21,S}
5    X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {7,S} {24,S}
6    X u0 p0 c0 s"fcc" m"terrace" {5,S} {10,S} {27,S}
7    X u0 p0 c0 s"hcp" m"terrace" {5,S} {11,S} {28,S}
8    X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {9,S} {31,S}
9    X u0 p0 c0 s"hcp" m"terrace" {8,S} {10,S} {34,S}
10   X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {9,S} {45,S}
11   X u0 p0 c0 s"bridge" m"terrace" {1,S} {7,S} {12,S} {49,S}
12   X u0 p0 c0 s"fcc" m"terrace" {11,S} {13,S} {50,S}
13   X u0 p0 c0 s"bridge" m"terrace" {1,S} {4,S} {12,S} {52,S}
14   X u0 p0 c0 s"ontop" m"terrace" {2,S} {15,S} {18,S} {20,S} {21,S} {23,S}
15   X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {17,S} {24,S}
16   X u0 p0 c0 s"fcc" m"terrace" {15,S} {20,S} {25,S}
17   X u0 p0 c0 s"hcp" m"terrace" {15,S} {23,S} {30,S}
18   X u0 p0 c0 s"bridge" m"terrace" {3,S} {14,S} {19,S} {31,S}
19   X u0 p0 c0 s"hcp" m"terrace" {18,S} {20,S} {32,S}
20   X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {19,S} {39,S}
21   X u0 p0 c0 s"bridge" m"terrace" {4,S} {14,S} {22,S} {52,S}
22   X u0 p0 c0 s"fcc" m"terrace" {21,S} {23,S} {53,S}
23   X u0 p0 c0 s"bridge" m"terrace" {14,S} {17,S} {22,S} {54,S}
24   X u0 p0 c0 s"ontop" m"terrace" {5,S} {15,S} {25,S} {27,S} {28,S} {30,S}
25   X u0 p0 c0 s"bridge" m"terrace" {16,S} {24,S} {26,S} {39,S}
26   X u0 p0 c0 s"hcp" m"terrace" {25,S} {27,S} {40,S}
27   X u0 p0 c0 s"bridge" m"terrace" {6,S} {24,S} {26,S} {45,S}
28   X u0 p0 c0 s"bridge" m"terrace" {7,S} {24,S} {29,S} {49,S}
29   X u0 p0 c0 s"fcc" m"terrace" {28,S} {30,S} {51,S}
30   X u0 p0 c0 s"bridge" m"terrace" {17,S} {24,S} {29,S} {54,S}
31   X u0 p0 c0 s"ontop" m"terrace" {8,S} {18,S} {32,S} {34,S} {36,S} {38,S}
32   X u0 p0 c0 s"bridge" m"terrace" {19,S} {31,S} {33,S} {39,S}
33   X u0 p0 c0 s"fcc" m"terrace" {32,S} {36,S} {42,S}
34   X u0 p0 c0 s"bridge" m"terrace" {9,S} {31,S} {35,S} {45,S}
35   X u0 p0 c0 s"fcc" m"terrace" {34,S} {38,S} {48,S}
36   X u0 p0 c0 s"bridge" m"terrace" {31,S} {33,S} {37,S} {49,S}
37   X u0 p0 c0 s"hcp" m"terrace" {36,S} {38,S} {51,S} 
38   X u0 p0 c0 s"bridge" m"terrace" {31,S} {35,S} {37,S} {54,S}
39   X u0 p0 c0 s"ontop" m"terrace" {20,S} {25,S} {32,S} {40,S} {42,S} {44,S}
40   X u0 p0 c0 s"bridge" m"terrace" {26,S} {39,S} {41,S} {45,S}
41   X u0 p0 c0 s"fcc" m"terrace" {40,S} {44,S} {46,S}
42   X u0 p0 c0 s"bridge" m"terrace" {33,S} {39,S} {43,S} {49,S}
43   X u0 p0 c0 s"hcp" m"terrace" {42,S} {44,S} {50,S}
44   X u0 p0 c0 s"bridge" m"terrace" {39,S} {41,S} {43,S} {52,S}
45   X u0 p0 c0 s"ontop" m"terrace" {10,S} {27,S} {34,S} {40,S} {46,S} {48,S}
46   X u0 p0 c0 s"bridge" m"terrace" {41,S} {45,S} {47,S} {52,S}
47   X u0 p0 c0 s"hcp" m"terrace" {46,S} {48,S} {53,S}
48   X u0 p0 c0 s"bridge" m"terrace" {35,S} {45,S} {47,S} {54,S}
49   X u0 p0 c0 s"ontop" m"terrace" {11,S} {28,S} {36,S} {42,S} {50,S} {51,S}
50   X u0 p0 c0 s"bridge" m"terrace" {12,S} {43,S} {49,S} {52,S}
51   X u0 p0 c0 s"bridge" m"terrace" {29,S} {37,S} {49,S} {54,S}
52   X u0 p0 c0 s"ontop" m"terrace" {13,S} {21,S} {44,S} {46,S} {50,S} {53,S}
53   X u0 p0 c0 s"bridge" m"terrace" {22,S} {47,S} {52,S} {54,S}
54   X u0 p0 c0 s"ontop" m"terrace" {23,S} {30,S} {38,S} {48,S} {51,S} {53,S}
"""
        slab_mol = Molecule().from_adjacency_list(slab_struct)
        rmol = Molecule().from_adjacency_list("""
1    X u0 p0 c0 {3,D}
2  * X u0 p0 c0 {4,R}
3    O u0 p2 c0 {1,D} {4,R}
4  * H u0 p0 c0 {2,R} {3,R}
""",check_consistency=False)
        site_map = {rmol.atoms[0]:slab_mol.atoms[32], rmol.atoms[1]:slab_mol.atoms[36]}
    
        out = get_full_mol(rmol,slab_mol,site_map)
    
        self.assertTrue(out.is_isomorphic(test_out,save_order=True))
    
    def test_get_labeled_full_TS_mol(self):
        full_mol = Molecule().from_adjacency_list("""1    X u0 p0 c0 s"ontop" m"terrace" {2,S} {5,S} {8,S} {10,S} {11,S} {13,S}
2    X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {4,S} {14,S}
3    X u0 p0 c0 s"fcc" m"terrace" {2,S} {8,S} {18,S}
4    X u0 p0 c0 s"hcp" m"terrace" {2,S} {13,S} {21,S}
5    X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {7,S} {24,S}
6    X u0 p0 c0 s"fcc" m"terrace" {5,S} {10,S} {27,S}
7    X u0 p0 c0 s"hcp" m"terrace" {5,S} {11,S} {28,S}
8    X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {9,S} {31,S}
9    X u0 p0 c0 s"hcp" m"terrace" {8,S} {10,S} {34,S}
10   X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {9,S} {45,S}
11   X u0 p0 c0 s"bridge" m"terrace" {1,S} {7,S} {12,S} {49,S}
12   X u0 p0 c0 s"fcc" m"terrace" {11,S} {13,S} {50,S}
13   X u0 p0 c0 s"bridge" m"terrace" {1,S} {4,S} {12,S} {52,S}
14   X u0 p0 c0 s"ontop" m"terrace" {2,S} {15,S} {18,S} {20,S} {21,S} {23,S}
15   X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {17,S} {24,S}
16   X u0 p0 c0 s"fcc" m"terrace" {15,S} {20,S} {25,S}
17   X u0 p0 c0 s"hcp" m"terrace" {15,S} {23,S} {30,S}
18   X u0 p0 c0 s"bridge" m"terrace" {3,S} {14,S} {19,S} {31,S}
19   X u0 p0 c0 s"hcp" m"terrace" {18,S} {20,S} {32,S}
20   X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {19,S} {39,S}
21   X u0 p0 c0 s"bridge" m"terrace" {4,S} {14,S} {22,S} {52,S}
22   X u0 p0 c0 s"fcc" m"terrace" {21,S} {23,S} {53,S}
23   X u0 p0 c0 s"bridge" m"terrace" {14,S} {17,S} {22,S} {54,S}
24   X u0 p0 c0 s"ontop" m"terrace" {5,S} {15,S} {25,S} {27,S} {28,S} {30,S}
25   X u0 p0 c0 s"bridge" m"terrace" {16,S} {24,S} {26,S} {39,S}
26   X u0 p0 c0 s"hcp" m"terrace" {25,S} {27,S} {40,S}
27   X u0 p0 c0 s"bridge" m"terrace" {6,S} {24,S} {26,S} {45,S}
28   X u0 p0 c0 s"bridge" m"terrace" {7,S} {24,S} {29,S} {49,S}
29   X u0 p0 c0 s"fcc" m"terrace" {28,S} {30,S} {51,S}
30   X u0 p0 c0 s"bridge" m"terrace" {17,S} {24,S} {29,S} {54,S}
31   X u0 p0 c0 s"ontop" m"terrace" {8,S} {18,S} {32,S} {34,S} {36,S} {38,S}
32   X u0 p0 c0 s"bridge" m"terrace" {19,S} {31,S} {33,S} {39,S}
33   X u0 p0 c0 s"fcc" m"terrace" {32,S} {36,S} {42,S} {55,S}
34   X u0 p0 c0 s"bridge" m"terrace" {9,S} {31,S} {35,S} {45,S}
35   X u0 p0 c0 s"fcc" m"terrace" {34,S} {38,S} {48,S}
36   X u0 p0 c0 s"bridge" m"terrace" {31,S} {33,S} {37,S} {49,S}
37   X u0 p0 c0 s"hcp" m"terrace" {36,S} {38,S} {51,S}
38   X u0 p0 c0 s"bridge" m"terrace" {31,S} {35,S} {37,S} {54,S}
39   X u0 p0 c0 s"ontop" m"terrace" {20,S} {25,S} {32,S} {40,S} {42,S} {44,S}
40   X u0 p0 c0 s"bridge" m"terrace" {26,S} {39,S} {41,S} {45,S}
41   X u0 p0 c0 s"fcc" m"terrace" {40,S} {44,S} {46,S}
42   X u0 p0 c0 s"bridge" m"terrace" {33,S} {39,S} {43,S} {49,S}
43   X u0 p0 c0 s"hcp" m"terrace" {42,S} {44,S} {50,S}
44   X u0 p0 c0 s"bridge" m"terrace" {39,S} {41,S} {43,S} {52,S}
45   X u0 p0 c0 s"ontop" m"terrace" {10,S} {27,S} {34,S} {40,S} {46,S} {48,S}
46   X u0 p0 c0 s"bridge" m"terrace" {41,S} {45,S} {47,S} {52,S}
47   X u0 p0 c0 s"hcp" m"terrace" {46,S} {48,S} {53,S}
48   X u0 p0 c0 s"bridge" m"terrace" {35,S} {45,S} {47,S} {54,S} {57,S}
49   X u0 p0 c0 s"ontop" m"terrace" {11,S} {28,S} {36,S} {42,S} {50,S} {51,S}
50   X u0 p0 c0 s"bridge" m"terrace" {12,S} {43,S} {49,S} {52,S}
51   X u0 p0 c0 s"bridge" m"terrace" {29,S} {37,S} {49,S} {54,S}
52   X u0 p0 c0 s"ontop" m"terrace" {13,S} {21,S} {44,S} {46,S} {50,S} {53,S}
53   X u0 p0 c0 s"bridge" m"terrace" {22,S} {47,S} {52,S} {54,S}
54   X u0 p0 c0 s"ontop" m"terrace" {23,S} {30,S} {38,S} {48,S} {51,S} {53,S}
55   O u0 p2 c0 {33,S} {56,S}
56   H u0 p0 c0 {55,S}
57   H u0 p0 c0 {48,S}""")
        rmol = Molecule().from_adjacency_list("""
        1     X u0 p0 c0 {3,S}
        2  *2 X u0 p0 c0 
        3  *3 O u0 p2 c0 {1,S} {4,S}
        4  *4 H u0 p0 c0 {3,S}
        5  *5 X u0 p0 c0 {6,S}
        6  *6 H u0 p0 c0 {5,S}
        """)
        labeled_mols = get_labeled_full_TS_mol(rmol,full_mol)
        test_out = Molecule().from_adjacency_list("""1     X u0 p0 c0 s"ontop" m"terrace" {2,S} {5,S} {8,S} {10,S} {11,S} {13,S}
2     X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {4,S} {14,S}
3     X u0 p0 c0 s"fcc" m"terrace" {2,S} {8,S} {18,S}
4     X u0 p0 c0 s"hcp" m"terrace" {2,S} {13,S} {21,S}
5     X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {7,S} {24,S}
6     X u0 p0 c0 s"fcc" m"terrace" {5,S} {10,S} {27,S}
7     X u0 p0 c0 s"hcp" m"terrace" {5,S} {11,S} {28,S}
8     X u0 p0 c0 s"bridge" m"terrace" {1,S} {3,S} {9,S} {31,S}
9     X u0 p0 c0 s"hcp" m"terrace" {8,S} {10,S} {34,S}
10    X u0 p0 c0 s"bridge" m"terrace" {1,S} {6,S} {9,S} {45,S}
11    X u0 p0 c0 s"bridge" m"terrace" {1,S} {7,S} {12,S} {49,S}
12    X u0 p0 c0 s"fcc" m"terrace" {11,S} {13,S} {50,S}
13    X u0 p0 c0 s"bridge" m"terrace" {1,S} {4,S} {12,S} {52,S}
14    X u0 p0 c0 s"ontop" m"terrace" {2,S} {15,S} {18,S} {20,S} {21,S} {23,S}
15    X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {17,S} {24,S}
16    X u0 p0 c0 s"fcc" m"terrace" {15,S} {20,S} {25,S}
17    X u0 p0 c0 s"hcp" m"terrace" {15,S} {23,S} {30,S}
18    X u0 p0 c0 s"bridge" m"terrace" {3,S} {14,S} {19,S} {31,S}
19    X u0 p0 c0 s"hcp" m"terrace" {18,S} {20,S} {32,S}
20    X u0 p0 c0 s"bridge" m"terrace" {14,S} {16,S} {19,S} {39,S}
21    X u0 p0 c0 s"bridge" m"terrace" {4,S} {14,S} {22,S} {52,S}
22    X u0 p0 c0 s"fcc" m"terrace" {21,S} {23,S} {53,S}
23    X u0 p0 c0 s"bridge" m"terrace" {14,S} {17,S} {22,S} {54,S}
24    X u0 p0 c0 s"ontop" m"terrace" {5,S} {15,S} {25,S} {27,S} {28,S} {30,S}
25    X u0 p0 c0 s"bridge" m"terrace" {16,S} {24,S} {26,S} {39,S}
26    X u0 p0 c0 s"hcp" m"terrace" {25,S} {27,S} {40,S}
27    X u0 p0 c0 s"bridge" m"terrace" {6,S} {24,S} {26,S} {45,S}
28    X u0 p0 c0 s"bridge" m"terrace" {7,S} {24,S} {29,S} {49,S}
29    X u0 p0 c0 s"fcc" m"terrace" {28,S} {30,S} {51,S}
30    X u0 p0 c0 s"bridge" m"terrace" {17,S} {24,S} {29,S} {54,S}
31    X u0 p0 c0 s"ontop" m"terrace" {8,S} {18,S} {32,S} {34,S} {36,S} {38,S}
32    X u0 p0 c0 s"bridge" m"terrace" {19,S} {31,S} {33,S} {39,S}
33    X u0 p0 c0 s"fcc" m"terrace" {32,S} {36,S} {42,S} {55,S}
34    X u0 p0 c0 s"bridge" m"terrace" {9,S} {31,S} {35,S} {45,S}
35    X u0 p0 c0 s"fcc" m"terrace" {34,S} {38,S} {48,S}
36    X u0 p0 c0 s"bridge" m"terrace" {31,S} {33,S} {37,S} {49,S}
37    X u0 p0 c0 s"hcp" m"terrace" {36,S} {38,S} {51,S}
38    X u0 p0 c0 s"bridge" m"terrace" {31,S} {35,S} {37,S} {54,S}
39    X u0 p0 c0 s"ontop" m"terrace" {20,S} {25,S} {32,S} {40,S} {42,S} {44,S}
40    X u0 p0 c0 s"bridge" m"terrace" {26,S} {39,S} {41,S} {45,S}
41    X u0 p0 c0 s"fcc" m"terrace" {40,S} {44,S} {46,S}
42    X u0 p0 c0 s"bridge" m"terrace" {33,S} {39,S} {43,S} {49,S}
43    X u0 p0 c0 s"hcp" m"terrace" {42,S} {44,S} {50,S}
44    X u0 p0 c0 s"bridge" m"terrace" {39,S} {41,S} {43,S} {52,S}
45    X u0 p0 c0 s"ontop" m"terrace" {10,S} {27,S} {34,S} {40,S} {46,S} {48,S}
46    X u0 p0 c0 s"bridge" m"terrace" {41,S} {45,S} {47,S} {52,S}
47    X u0 p0 c0 s"hcp" m"terrace" {46,S} {48,S} {53,S}
48 *5 X u0 p0 c0 s"bridge" m"terrace" {35,S} {45,S} {47,S} {54,S} {57,S}
49    X u0 p0 c0 s"ontop" m"terrace" {11,S} {28,S} {36,S} {42,S} {50,S} {51,S}
50    X u0 p0 c0 s"bridge" m"terrace" {12,S} {43,S} {49,S} {52,S}
51    X u0 p0 c0 s"bridge" m"terrace" {29,S} {37,S} {49,S} {54,S}
52    X u0 p0 c0 s"ontop" m"terrace" {13,S} {21,S} {44,S} {46,S} {50,S} {53,S}
53    X u0 p0 c0 s"bridge" m"terrace" {22,S} {47,S} {52,S} {54,S}
54    X u0 p0 c0 s"ontop" m"terrace" {23,S} {30,S} {38,S} {48,S} {51,S} {53,S}
55 *3 O u0 p2 c0 {33,S} {56,S}
56 *4 H u0 p0 c0 {55,S}
57 *6 H u0 p0 c0 {48,S}""")
        self.assertTrue(test_out.is_isomorphic(labeled_mols[0],generate_initial_map=True,save_order=True))
                        
if __name__ == '__main__':
    unittest.main()