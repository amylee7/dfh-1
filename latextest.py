""" Unit Test for latex.py"""

from Latex import *
from tangle_ import TANGLE, StrandDiagram, Simple_Strand, Idempotent, StrandAlgebra, TangleModule
import unittest
from utility import F2

class LatexTest(unittest.TestCase):
    def testPrintDAStructure(self):
        
        dic = {(4,4):(3.5,4), (3.5,4):(3,4), (3,2):(3.5,2), (3.5,2):(4,3), (4,2):(3.5,3), (3.5,3):(3,3), (3,1):(3.5,1), (3.5,1):(4,1)}
#        dic = {(4,4):(3.5,4),(3.5,4):(3,4),(4,2):(3.5,3), (3.5,3):(3,3), (3,2):(3,1)}
        tang = TANGLE(dic)
        strands =  (((3,4),(2,1),(1,2)),((3,2),(0,1)))
        right_alg = StrandAlgebra(F2, tang, False)
        a1 = Simple_Strand(right_alg, False, ((2,2),(1,4)))
        CT = TangleModule(F2, tang) #CT(Ti)
        sd = StrandDiagram(CT,strands) # strand diagram
        f = open("latex_output.txt", "w")
        f.write(beginDoc())
        f.write(showTensor(sd,a1))
        f.write(endDoc())
        f.close()

if __name__ == "__main__":
    unittest.main()
    

