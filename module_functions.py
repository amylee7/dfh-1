from tangle_ import Simple_Strand, StrandAlgebra,TANGLE,Strands, Idempotent, StrandDiagram
from utility import mod_helper, mod_between
import numpy as np

import operator
#from algebra import DGAlgebra, Element, Generator, Tensor, TensorGenerator
import statistics as st
from utility import orientation,orientation_i, complement,generate_subset,\
 doescross, doescross_simple, intersections, simple_intersections, F2
from Algebra import DGAlgebra, Element, Generator, Tensor, TensorGenerator
from Algebra import E0
from utility import get_domain, get_range, generate_bijections, combinations, \
                get_domain_dict, get_range_dict, doescross_simple_rc,in_between_list, \
                reorganize_sign, reorganize_sign_2,dict_shift, get_start_dict, get_end_dict,\
                dict_shift_double, doescross_bool, mod_between, mod_helper

def del_l_raw(gen):
    ''' del_L map that sends CT-(T_i) to  A-(dRT) * I(-dRT) * CT-(T_i)
    gen is an element of CT-(T_i).'''
    if not isinstance(gen, StrandDiagram):
        return NotImplemented
    prod_raw = d_m_raw_A(gen)
    if prod_raw is None:
        return E0
    if self.ring is F2: # cb and change
        return prod_raw.elt()   
    else:
        pass

def m2(gen, alg_gen):
    ''' m2 map that sends CT-(T) * I(-dRT) * A-(dRT) to CT-(T).
    * gen is a generator in CT-(Ti)
    * alg_gen is a generator of  A-(dRTi) 
    gen1 here refers to the strand diagram generator of right most elementary tangle T_i '''  
    if not isinstance(gen, StrandDiagram):
        return NotImplemented
    if not isinstance(alg_gen, Simple_Strand):
        return NotImplemented
    assert gen.tangle == alg_gen.p_tangle, " Not compatible. Different Tangle."
    assert alg_gen.is_left == False, " Not compatible, Algebra takes -dLT."    
    if not gen.strands.rightCompatible(gen, alg_gen):
        return E0
    prod_raw = m2_raw(gen, alg_gen)
    if prod_raw is None:
        return E0
    if self.ring is F2: # cb and change
        return prod_raw.elt()
    else:
        pass

def m2_raw(gen, alg_gen):
    ''' If gen and alg_gen can be multiplied (i.e.,the strands pairs matching up), return the generator in CT-(T) that is
    their product. Otherwise return None. gen1 and alg_gen 2 are compatible by default.''' 
    if mod_6_m2(gen, alg_gen):
        return None
    else: #concantenate
        l_raw = gen.strands.data[1]
        r_raw = alg_gen.strands.data[0]
        new = conc_strands(l_raw, r_raw)
        return gen.replace_2([[],list(l_raw)],[[],new])
    
def mod_6_m2(gen, alg_gen):
    ''' Mod 6 relations for m2'''    
    l_tangle = [gen.tangle.orient_right_rhalf, gen.tangle.orient_left_rhalf]
    r_tangle = [alg_gen.tangle.orient_right_lhalf, alg_gen.tangle.orient_left_lhalf]  
    l_strands = gen.strands.right_converted
    r_strands = alg_gen.strands.left_converted
    c_l = gen.strands.strandCrossing_coord(False)
    c_r = gen.strands.strandCrossing_coord(True)
    return len(mult_two_halfs(l_tangle, l_strands, r_tangle, r_strands)) > 0 \
                or len(cross_twice(l_strands, r_strands, c_l, c_r)) > 0 
    
def d_plus_raw_B(gen):
    ''' d_+ map that smooths a black-black crossing of a generator `gen`
    of CT-(T_i)) contained in right halfs ( i- 1/2, i)
    Returns a list of elements of the form ((s1, s2), d_plus_raw_term), where
    s1 < s2 are pairs of strands that d_plus is applied, and
    diff_term is a generator in d_plus() obtained by uncrossing these two
    strands. Together they specify all terms in gen.d_plus(). '''    
    # applying it to line Bj
    lst = []
    l_strands = gen.strands.left_converted 
    r_strands = gen.strands.right_converted
    types = mod_helper(l_strands, r_strands)
    t_r = gen.tangle.r_pairs_wc
    r_crossings = types[4] # type 4 : right crossings
    for s1, s2 in r_crossings: 
        is_crossed = False
        # check double black - black crossing
        for strands in r_strands: 
            if mod_between(strands, (s1[0],s2[0]),(s1[1],s2[1]), False) == 1:
                is_crossed = True
                break
        for tangle in t_r:# tangle double cross a strand
            if mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]), False) ==1:
                is_crossed = True
                break
        if not is_crossed:
                lst.append(mod_helper_2(gen, s1, s2, is_crossed))
    return lst

def d_minus_raw(gen, is_B):
    '''d_- map that introduces a black-black crossing to a generator 'gen`
    of CT-(T_i) contained in left halfs, (i, i+1/2).
    if is_B == True, apply it to B_j
    if is_B == False apply it to A_j 
    
    Returns a list of elements of the form ((s1, s2), d_minus_raw_term), where
    s1 < s2 are pairs of strands that d_minus is applied, and
        diff_term is a generator in d_minus() obtained by crossing these two
        strands. Together they specify all terms in gen.d_minus(). ''' 
    lst = [] # applying it to line Bj
    l_strands = gen.strands.left_converted 
    r_strands = gen.strands.right_converted
    types = mod_helper(l_strands, r_strands)
    t_l = gen.tangle.l_pairs_wc
    # type 1 : left horizontals
    horizontals = types[1]
    for s1, s2 in horizontals:
        is_crossed = False
        for strands in r_strands: # check double black - black crossing
            if mod_between(strands, (s1[0],s2[0]),(s1[1],s2[1]),True)==0:
                is_crossed = True  
                break
        for tangle in t_l:# tangle double cross a strand
            if mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),True) == 0:
                is_crossed = True
                break
        if not is_crossed:
                lst.append(mod_helper_2(gen, s1, s2, is_crossed))
    return lst

def d_m_raw_A(gen, algebra, alg_left):
    '''dm map that picks two pair of points along  A_j, bewteen 
    Exchanges two ends of the corresponding pair of black strands.
    if `alg_left` == True, A(-dLT_i) * CT(T_ii)
    Returns a list of elements of the form ((s1, s2), d_m_raw_term), where
    s1 < s2 are pairs of strands that d_m is applied, and
        diff_term is a generator in d_plus() obtained by uncrossing or  crossing
        these two strands. Together they specify all terms in d_m()'''
    
     # So far only implemented for A(-dLT_i) * CT(T_ii)
    
    lst =  []  # CB
    l_strands = []
    r_strands = []
    t_l = []
    t_r = []
    types = mod_helper(l_strands, r_strands)
    # left no cross strands 
    l_horizontals = types[1]
    for s1, s2 in l_horizontals:
        is_crossed = False
        for tangle in t_l: # left half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]), True)) == 1:
                is_crossed = True
                break
                break
        for tangle in t_r: # right half tangles
            x = (tuple(np.add(s1[0],(1,0))), tuple(np.add(s2[0],(1,0))))
            if mod_between(tangle,x, False ) != -2:
                is_crossed = True
                break
        if not is_crossed:
            pass
    # right cross strands
    r_cross = types[4]
    for s1, s2 in r_cross:
        is_crossed = False
        for tangle in t_l: # left half tangles
            x = (tuple(np.subtract(s1[1],(1,0))), tuple(np.subtract(s2[1],(1,0))))
            if mod_between(x,(s1[0],s2[0]),True) != -2:
                is_crossed = True
                break
        for tangle in t_r: # right half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),False)) == 1:
                is_crossed = True
                break
    # left below right
    l_below_r = types[5]
    for s1, s2 in l_below_r:
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: #left half tangles
            if mod_between(tangle, left_pair, mid_pair, True) == -1:
                is_crossed = True
                break
        for tangle in t_r: # right half tangles
            if -2 < mod_between(tangle, mid_pair, right_pair, False) < 1:
                is_crossed = True
                break       
    # left above right strands
    l_above_r = types[6]
    for s1, s2 in l_above_r:
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: # left half tangle
            if mod_between(tangle, left_pair, mid_pair, True) == 1:
                is_crossed = True
                break
        for tangle in t_r: # right half tangle
            if mod_between(tangle, mid_pair, right_pair, False) > -1:
                is_crossed = True
                break
    return lst

def d_m_raw_B(gen):
    ''' dm map that picks two pair of points along Bj.
    Exchanges two ends of the corresponding pair of black strands of generator
    `gen` in CT(Ti)
        Returns a list of elements of the form ((s1, s2), d_m_raw_term), where
    s1 < s2 are pairs of strands that d_m is applied, and
        diff_term is a generator in d_plus() obtained by uncrossing or  crossing
        these two strands. Together they specify all terms in d_m(). '''
    lst = []
    l_strands = gen.strands.left_converted 
    r_strands = gen.strands.right_converted
    types = mod_helper(l_strands, r_strands)
    t_l = gen.tangle.l_pairs_wc 
    t_r = gen.tangle.r_pairs_wc 
    # right no cross strands
    r_horizontals = types[3]
    for s1,s2 in r_horizontals:
        is_crossed = False
        for tangle in t_r: # right half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),False)) == 1:
                is_crossed = True
                break
        for tangle in t_l: #left half tangles
            x = (tuple(np.subtract(s1[1],(1,0))), tuple(np.subtract(s2[1],(1,0))))
            if mod_between(tangle, x,(s1[0],s2[0]), True) != -2:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2, is_crossed))
    # left cross strands
    l_cross = types[2]
    for s1,s2 in l_cross:
        is_crossed = False
        for tangle in t_l: #left half tangles
            if abs(mod_between(tangle, (s1[0],s2[0]),(s1[1],s2[1]),True)) == 1:
                is_crossed = True
                break
        for tangle in t_r: # right half tangles
            x = (tuple(np.add(s1[0],(1,0))), tuple(np.add(s2[0],(1,0))))
            if mod_between(tangle,(s1[1],s2[1]),x, False) != -2:
                is_crossed = True 
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2, is_crossed))
    # left above right strands
    l_above_r = types[6]
    for s1, s2 in l_above_r:
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: # left half tangles
            if -2 < mod_between(tangle, left_pair, mid_pair, True) < 1:
                is_crossed = True
                break
        for tangle in t_r: # right_half tangles
            if mod_between(tangle, mid_pair, right_pair, False) == -1:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2, is_crossed)) 
    # left below  right strands
    l_below_r = types[5]      
    for s1,s2 in l_below_r: #s1 on the left, s2 on the right
        is_crossed = False
        left_pair = (tuple(np.subtract(s2[0],(0.5,0))), s1[0])
        mid_pair = (s2[0],s1[1])
        right_pair = (tuple(np.add(s1[1],(0.5,0))),s2[1])
        for tangle in t_l: # left half tangles
            if mod_between(tangle, left_pair, mid_pair, True) > -1 :
                is_crossed = True
                break
        for tangle in t_r:
            if mod_between(tangle, mid_pair, right_pair, False) == 0:
                is_crossed = True
                break
        if not is_crossed:
            lst.append(mod_helper_2(gen, s1, s2, is_crossed))
    return lst

def mod_helper_2(gen,s1, s2, is_crossed):
    '''returns the new generator if is_crossed is false '''
    if not is_crossed:
        if isinstance(s1[0][0], int):
            sd_1 = gen.strands.get_strand_index(s1, True)
            sd1_left = True
        else:
            sd_1 = gen.strands.get_strand_index(s1, False)
            sd1_left = False  
        if isinstance(s2[0][0], int):
            sd_2 = gen.strands.get_strand_index(s2, True)
            sd2_left = True
        else:
            sd_2 = gen.strands.get_strand_index(s2, False)
            sd2_left = False    
        if sd1_left and sd2_left:
            new_strand = ((sd_1[0],sd_2[1]),(sd_2[0],sd_1[1]))
            return (sd_1, sd_2), gen.replace((sd_1, sd_2), new_strand, True)
        elif not sd1_left and not sd2_left:
            new_strand = ((sd_1[0],sd_2[1]),(sd_2[0],sd_1[1]))
            return ((sd_1, sd_2), gen.replace((sd_1, sd_2), new_strand, False))
        elif sd1_left and not sd2_left:
            new_strand = ((sd_1[0],sd_2[0]),(sd_1[1],sd_2[1]))
            return ((sd_1, sd_2), gen.replace_2((sd_1, sd_2), new_strand))
        else:
            new_strand = ((sd_1[0],sd_2[0]),(sd_1[1],sd_2[1]))
            return ((sd_1, sd_2, gen.replace_2((sd_1, sd_2), new_strand)))
    else:
        return None

def helper(tang_left, tang_right, pair_1, pair_2, pair_3, is_left, option):
    ''' Helper methods for diff methods. Assumes compatibility. Given a generator x, taking the is_left side, 
    exchanges pairs = (a,b), it seems whether the 
    whether it violates `option`. Refer to options.jpg
    for options.  
    * tang_left is dictionary object of pairs 
    * tang_right is dictionary object of pairs 
    * Pair_1 is the coordinate of two pairs of strand on the left
    * Pair_2 is the coordinate of two pairs of strand in the middle
    * Pair_3 is the coordinate of two pairs of strand in the middle.'''
    mod = False
    if is_left: # only check mod relations 1,2,3
        for k,v in tang_left.items():
            if option ==1:
                if mod_helper((k,v), pair_1, pair_2, is_left) == 1: 
                    mod = True
            elif option == 2:
                if mod_helper((k,v), pair_1, pair_2, is_left) == 0:
                    mod = True
            elif option == 3:
                if mod_helper((k,v), pair_1, pair_2, is_left) == -1: 
                    mod = True
            else:
                raise AssertionError("Something is wrong -- options arent either 1, 2 or 3 ")      
    else: # only check mod relations 4,5,6
        for k,v in tang_right.items():
            if option == 4:
                if mod_helper((k,v), pair_2, pair_3, is_left) == 1:
                    mod = True
            elif option ==5:
                if mod_helper((k,v), pair_2, pair_3, is_left) == 0: 
                    mod = True
            elif option == 6:
                if mod_helper((k,v), pair_2, pair_3, is_left) == -1: 
                    mod = True    
            else:
                raise AssertionError("Something is wrong about `is_left`-- options arent either 4, 5, 6" )
    return mod

## helper methods for mod 6() 
    def mult_two_halfs(self,left_tangle, left_strands, right_tangle, right_strands):
        ''' Multiplies two half tangles T1 and T2, and mod out first two
        relations in PetKova paper Figure 6. 
        left_tangle : [dict1, dict2] right half of tangle T1, array of dictionary objects of coordinates
                      dict1 : orienting left, dict2: orienting right
        right_tangle: [dict1, dict2] left half of tangle T2, a dictionary object of coordinates
        left_strand: right half of tangle T1 strands, a dictionary object of coordinates
        right_strand: left half of tangle T1 strands, a dictionary object of coordinates
        
        pairs is an array object, which elements of form [start tangle, start strand]
        used for algebra multiplication and for other things later( d+, dm)'''
        pairs = []
        # tangle orienting right      
        t1 = left_tangle[0]
        t2 = right_tangle[0]   
        if(self.is_left): # left side of the tangle
            left_strands = dict_shift(left_strands,True)# shift left_strand to the right
        else:
            right_strands = dict_shift(right_strands,False) # shift left_strand to the right      
        for k,v in t1.items():
            for a,b in left_strands.items():
                if doescross((k,v),(a,b)):
                    new_tangle = (v,t2[v]) # on the right
                    new_strand = (b, right_strands[b]) # on the right
                    if doescross(new_tangle, new_strand):
                        pairs.append([k,a])       
        # tangle orienting left
        t1 = left_tangle[1]
        t2 = right_tangle[1]               
        for k,v in t2.items():
            for b,c in right_strands.items():
                if doescross((v,k),(b,c)):
                    new_tangle = (t1[v], v) # on the left        
                    # search for a strand in left_strand, that ends in b
                    for key,value in left_strands.items():
                        if value == b:
                            new_strand = (key,b) # get the left coordinate     
                    if doescross(new_tangle, new_strand):
                        pairs.append([k, c])
        return pairs

    def cross_twice(self,left_strands, right_strands, c_l, c_r): # cb and change it to true or false later if uncessary
        '''Given left_strand, and right_strand, returns double crossings,
        assumes the strands are multiplicable. 
        c_1 is strandCrossing(False) called from the left_strands
        c_2 is the strandCrossing(True) called from right_strands'''
        # check if compatible cb and remove
        if(self.is_left): # left side of  the tangle
            left_strands = dict_shift(left_strands,True)# shift left_strand to the right
            c_l = dict_shift_double(c_l ,True)
        else:
            right_strands = dict_shift(right_strands,False) # shift left_strand to the right 
            c_r = dict_shift_double(c_r ,False)        
        if not (get_range_dict(left_strands) == get_domain_dict(right_strands)):
            raise TypeError ("strand boundaries don't line up")         
        double_crossings = []
        for k in list(c_l.items()):
            for v in list(c_r.items()):
                if get_end_dict(k) == get_start_dict(v):
                    double_crossings.append(k)
        return double_crossings
    
#### TEST CODE ####
p1 = ((1,0.5),(1,4.5))
p2 = ((1.5,1.5), (1.5, 3.5))
p3 = ((2,1.5),(2,4.5))
t_left = {(1,7):(1.5,3), (1,2):(1.5,3)}
t_right = {(1.5,3):(2,5)}

print(helper(t_left,t_right, p1, p2, p3, False, 6 ))
         