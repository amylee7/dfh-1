#Bohua Zhan @https://github.com/bzhan/bfh_python

# Adapted from Bzhan Code are as follows: 
#NamedObject, SummableDict, fracToInt, memorize, memorizeHash, safeMultiply, F2

from fractions import gcd
from numbers import Number
import itertools as it
from numpy import *
import numpy as np

# adapted from bzhan

class Ring:
    def convert(self, data):
        """Try to convert data to an element of this ring."""
        raise NotImplementedError("convert function not specified for ring.")

class RingElement(Number):
    pass


class ModNRing(Ring):
    """The ring Z/nZ."""
    def __init__(self, n):
        self.n = n
        self.zero = self.convert(0)
        self.one = self.convert(1)

    def add(self, elt1, elt2):
        elt1, elt2 = self.convert(elt1), self.convert(elt2)
        if elt1 is NotImplemented or elt2 is NotImplemented:
            return NotImplemented
        return ModNElement(self, (elt1.val+elt2.val)%self.n)

    def multiply(self, elt1, elt2):
        elt1, elt2 = self.convert(elt1), self.convert(elt2)
        if elt1 is NotImplemented or elt2 is NotImplemented:
            return NotImplemented
        return ModNElement(self, (elt1.val*elt2.val)%self.n)

    def __eq__(self, other):
        return self.n == other.n

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.n, "ModNRing"))

    def convert(self, data):
        """Try to convert data to an element of this ring."""
        if isinstance(data, ModNElement) and data.parent == self:
            return data
        if isinstance(data, int):
            return ModNElement(self, data % self.n)
        return NotImplemented

class ModNElement(RingElement):
    """An element in a ring Z/nZ."""
    def __init__(self, parent, val):
        self.parent = parent
        self.val = val

    def __str__(self):
        return str(self.val)

    def __repr__(self):
        return str(self.val)

    def __add__(self, other):
        return self.parent.add(self, other)

    def __radd__(self, other):
        return self.parent.add(self, other)

    def __mul__(self, other):
        return self.parent.multiply(self, other)

    def __rmul__(self, other):
        return self.parent.multiply(self, other)

    def __eq__(self, other):
        """Can compare to integer 0 or 1."""
        if isinstance(other, int) and (other == 0 or other == 1):
            return self.val == other
        else:
            return self.val == other.val

    def invertible(self):
        """Returns whether this element is invertible in the ring."""
        return gcd(self.val, self.parent.n) == 1

    def inverse(self):
        """Returns the inverse of this element. Must be invertible"""
        # Currently only implemented for n = 2
        assert self.parent.n == 2
        return self

class Integer(Ring):
    """The ring Z."""
    def convert(self, data):
        """Try to convert data to an element of this ring."""
        assert isinstance(data, int)
        return data

class IntegerElement(RingElement, int):
    """An element in a ring Z."""
    def __new__(cls, parent, val):
        return int.__new__(cls, val)

    def __init__(self, parent, val):
        self.parent = parent


F2 = ModNRing(2)
ZZ = Integer()

# Constants for positive and negative orientation
POS, NEG = 1, -1

class NamedObject:
    """Provides functionality for an object to be described by name. If this is
    listed as a parent class, an object will use name in equality comparisons,
    hash functions, and string outputs.
    """
    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return not (self == other)

    #memorizeHash
    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return str(self.name)

    def __repr__(self):
        return str(self.name)


#an analogous version Free Modules by extending dictionary data structure
class SummableDict(dict):
    ''' An extension of Dictionary class that supports sum and multiplication:
        Works like free module with generators as keys, and coefficients as values. 
        Zeros are automatically thrown away'''
        
    def __add__(self, other):
        ''' does not modify the invoking object '''
        return _dictAddTo(type(self)(), [self, other])
    
    def __iadd__(self, other):
        '''typical add operation for dictionaries: add other onto self
              modifies the invoking object.'''
        return _dictAddTo(self, other)
    
    def __sub__(self, other):
        
        # ORIGINAL : return _dictAddTo(type(self)(), [self, -1*other])
        return _dictAddTo(type(self)(), [self, _dictMult(other, -1)])
    
    def __isub__(self, other):
        '''typical subtract operation for dictionaries: subtract other from self'''
        #ORIGINAL return _dictAddTo(self, -1*other)
        return _dictAddTo(self, _dictMult(other, -1))
    
    def accumulate(self, lst):
        """Similar to +=, except returns the sum."""
        return _dictAddTo(self, lst)
       
    def __mul__(self, other):
        '''multiply by scalar'''
        return _dictMult(self, other)
    
    def __rmul__(self, other): #useless?
        # dk what the difference is against _mult # cb
        return _dictMult(self, other)
    
    def __eq__(self, other):
        ''' Comparison to 0 is a test for empty dictionary.'''
        if isinstance(other, int) and other == 0:
            return len(self) == 0
        else:
            return dict.__eq__(self, other)
    
    def __ne__(self, other): 
        '''returns whether self is equal to other'''
        return not (self == other)
    
    def copy(self): #useless?
        '''Copy function that preserves type.'''
        return type(self)(self)
        
    def translateKey(self, key_map):
        '''Translate keys of self using key_map. 
        keys in self must appear in key_map and key named will be changed to the 
        corresponding values (i.e the new key names) in key_map dictionary'''
        
        return type(self)([(key_map[k], v) for k, v in self.items()])
    
    def getElt(self):
        """Returns an arbitrary key from this dictionary. Must be non-empty."""
        #return iter(self).next() # in python 2
        return iter(self).__next__() # in python 3 

#dictionary helper methods:
    
def tolist(obj):
    """Force obj into a list."""
    if isinstance(obj, list):
        return obj
    else:
        return [obj]    
#adding dictionaries

def _dictAddTo(dict1, dict2):
    """Add dict2 onto dict1 in place. If dict2 is a list of dictionaries,
    add each element of
    dict2 on dict1 in place.
    """
    dict2 = [curdict.copy() for curdict in tolist(dict2) if curdict != 0] # removes all empty dictionary
    if dict1 == 0:
        if len(dict2) == 0:
            return dict1
        else:
            dict1 = dict2[0]
            dict2 = dict2[1:]
    for curdict in dict2:
        assert type(dict1) == type(curdict), "Incompatible types: %s, %s" % \
            (str(type(dict1)), str(type(curdict)))
        for k, v in curdict.items():
            #if dict1.has_key(k): 
                # has key is removed in python 3
            if k in dict1:
                dict1[k] += v
                if dict1[k] == 0:
                    del dict1[k]
            else:
                dict1[k] = v
    return dict1
    
#multiplying dictionaries
def _dictMult(dict1, scalar):
    '''Return a new dictionary with same type as self, the same keyes, and each value mlutiplied by scalar'''
    
    if not isinstance(scalar, Number):
        return NotImplemented
    
    result = type(dict1)((k, scalar * v) for k, v in dict1.items() if
                             scalar * v != 0)
    return result

def safeMultiply(a, b):
    '''multiply two sides using __mul__ and __rmul__. 
    Return NotImplemented if both fails.'''
    
    try:
        prod = a.__mul__(b)
    except TypeError:
        prod = NotImplemented
        
    if prod is NotImplemented:
        try:
            prod = b.__rmul__(a)
        except TypeError:
            prod = NotImplemented
    return prod

# ------- practice ------- dict codes

def __add__(self, other):
    '''typical add operation for dictionaries: add other onto self'''
    print("THIS IS WHAT IT MEANS")
    a = type(self)()
    print(a)
    print(type(a))
    return _dictAddTo(type(self)(), [self, other])

def __iadd__(self, other):
    '''like __add__ but other can be a list of dictionaries.'''
    return _dictAddTo(self, other)

def __sub__(self, other):
    '''typical subtract operation for dictionaries: subtract other from self'''
    return _dictAddTo(type(self)(), [self, _dictMult(other, -1)])#diff

def __isub__(self, other):
    '''like __sub__ but other can be a list of dictionaries.'''
    return _dictAddTo(self, _dictMult(other, -1)) #? diff

def accumulate(self, lst):
    """Similar to +=, except returns the sum."""
    return _dictAddTo(self, lst)
   
def __mul__(self, other):
    '''multiply by scalar''' #?
    return _dictMult(self, other)

def __rmul__(self, other): #useless?
    # dk what the difference is against _mult #?
    return _dictMult(self, other)

def __ne__(self, other):
    '''returns whether self is equal to other'''
    return not (self == other)

def copy(self): #useless?
    '''Copy function that preserves type.'''
    return type(self)(self)
    
def translateKey(self, key_map):
    '''Translate keys of self using key_map. 
    keys in self must appear in key_map and key named will be changed to the 
    corresponding values (i.e the new key names) in key_map dictionary'''
    return type(self)([(key_map[k], v) for k, v in self.items()])

def getElt(self):
    """Returns an arbitrary key from this dictionary. Must be non-empty."""
    #return iter(self).next() # in python 2
    return iter(self).__next__() # in python 3 

def complement(tup1, tup2):
    ''' Given ordered tup1 = [n] and returns the tup1 \ tup2,
    given that tup2 is contained in tup1''' 

    return tuple(set(tup1).difference(set(tup2)))


def orientation(pair, left_half = True):
    ''' Given a pair of coordinates ((x_1,y_1), (x_2,y_2)), returns 
    list {1} if from left to right 
    list {-1} if  right to left
    (1, -1) or (1,1) depending on cup or cap '''
    
    # left_half is `True` if is left half of the tangle i to i+1/2
    # left_half is `False` if it is right half of the tangle i+1/2 to i+1

    x_1,y_1 = pair[0]
    x_2,y_2 = pair[1]
    
    if x_1 < x_2: # going left to right
        return (1)
    elif x_1 > x_2:
        return (-1)

    else: # cups or caps 
        if left_half == True:
            if y_1 > y_2:
                return (-1, 1)
            else:
                return (1, -1)
        else:#right half
            if y_1 > y_2:
                return (1, -1)
            else:
                return (-1, 1)

            
def orientation_i(pair, left_half = True):
     ''' Improved version of `orientation` method above, but actually returns 
     a dictionary object {(x1,y1): +1, (x2, y2): +- 1,'''
     x_1,y_1 = pair[0]
     x_2,y_2 = pair[1]
     pair_0 = (x_1,y_1)
     pair_1 = (x_2,y_2)

    # IF LEFT HALF
    
     if left_half == True:
        if x_1 < x_2: # going left to right
             return {pair_0 :1} # cb
        elif x_1 > x_2:
            return {pair_1 :-1}
        else: # caps
            if y_1 > y_2:
                return {pair_1: -1, pair_0: 1}
            else:
                return {pair_0: 1, pair_1: -1} 
    # IF RIGHT HALF
     else:
        assert left_half == False # assert error if not right half
        if x_1 < x_2: # going left to right
             return {pair_1 :1} # cb
        elif x_1 > x_2:
            return {pair_0 :-1}
        else : #caps
            if y_1 > y_2:
                 return {pair_1: 1, pair_0: -1}
            else:
                 return {pair_0: -1, pair_1: 1}    


def doescross_simple(left,right):
    '''given a pair of points left = (a_1, b_1),right = (a_2, b_2)
    it does tells us whether they intersect''' 
    
    a = ((1, left[0]),(2,left[1]))
    b = ((1, right[0]),(2, right[1]))
    
    def line(p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0]*p2[1] - p2[0]*p1[1])
        return A, B, -C
    
    def intersection(L1, L2):
        D  = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        if D != 0:
            x = Dx / D
            y = Dy / D
            return x,y
        else:
            return False
        
    L1 = line(a[0],a[1])
    L2 = line(b[0],b[1])
    
    R = intersection(L1, L2)
    #R = seg_intersect(a,b)
    if R == False:
        R = None
    elif R[0] < 1 or R[0] > 2: 
        R = None
    elif not(in_between(left[0], left[1],R[1]) and in_between(right[0], right[1], R[1])): # Check Y coordinates
        R = None
    else:
        pass
    
    if R != None:
        return True
    else:
        return False

def doescross_simple_rc(left,right):
    '''given a pair of points left = (a_1, b_1),right = (a_2, b_2)
    it does tells us whether they intersect and RETURNS coordinates''' 
    
    a = ((1, left[0]),(2,left[1]))
    b = ((1, right[0]),(2, right[1]))
    
    def line(p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0]*p2[1] - p2[0]*p1[1])
        return A, B, -C
    
    def intersection(L1, L2):
        D  = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        if D != 0:
            x = Dx / D
            y = Dy / D
            return x,y
        else:
            return False
        
    L1 = line(a[0],a[1])
    L2 = line(b[0],b[1])
    
    R = intersection(L1, L2)
    #R = seg_intersect(a,b)
    if R == False:
        R = None
    elif R[0] < 1 or R[0] > 2: 
        R = None
    elif not(in_between(left[0], left[1],R[1]) and in_between(right[0], right[1], R[1])): # Check Y coordinates
        R = None
    else:
        pass
    
    return R

def in_between(a,b, x):
    ''' returns True if x is in between a and b. 
        a > b or b < a may happen. Error returned if a = b'''
        
    if a==b and a != x:
        raise ValueError("a = b. Something is wrong.")
    if a == b and a == x:
        return True
    elif a>b:
        return b < x and x < a
    else:
        return a < x and x < b
    
def in_between_list(arr, x):
    ''' Given an array `arr`, returns true if there exists two values a1, a2 
    in `arr` such that x is in between arr
    if arr size is one, returns type error'''
    if not isinstance(arr,list):
        return TypeError("`arr` is not an array. Wrong type of input")
    elif len(arr) == 1:
        return TypeError("`arr` is of size 1. cannot use this method.")
    else:
        arr.sort()
        maximum = max(arr)
        minimum = min(arr)
        if x > maximum or x < minimum:
            return False
        else:
            return True
        
def reorganize_sign(cross):
    '''returns [a1, a2, b1, b2]
    Reorient points such that a1 is top left, then b1, a2, b2 in counter-
    clockwise orientation'''
    
    # ex if cross = ((4,2),(1,3))
    s1 = cross[0][0] #4
    s2 = cross[1][0] #1
    
    if s1 > s2:
        a1 = cross[0][0]
        a2 = cross[0][1]
        b1 = cross[1][0]
        b2 = cross[1][1]
    
    else: # s1 < s2
        a1 = cross[1][0]
        a2 = cross[1][1]
        b1 = cross[0][0]
        b2 = cross[0][1]
    
    return [a1,a2,b1,b2]
def doescross(a,b):
    '''given a pair of lines = ((x_1, y_1), (x_2, y_2)), ((x_3, y_3), (x_4, y_4))
    goes the intersection point, if one exists between, if not return None''' 
    
    def line(p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0]*p2[1] - p2[0]*p1[1])
        return A, B, -C
    
    def intersection(L1, L2):
        D  = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        if D != 0:
            x = Dx / D
            y = Dy / D
            return x,y
        else:
            return False
        
    L1 = line(a[0],a[1])
    L2 = line(b[0],b[1])
    
    R = intersection(L1, L2)

    
    x_1, y_1 = a[0]
    x_2, y_2 = a[1]
    x_3, y_3 = b[0]
    x_4, y_4 = b[1]
    

    if R == False:
        R = None
    elif not(in_between(x_1, x_2, R[0]) and in_between(x_3, x_4, R[0]) and \
           in_between(y_1, y_2, R[1]) and in_between(y_3, y_4, R[1])):
        R = None
    else:
        pass

    return R


def combinations(n, m):
    ''' retyrbs an array of all possible combinations of n and m
    for example, if n = 1 , n = 2, returns, if WLOG n < m
    [ (0,0),(0,1),.....(n, m-1 ), (n,m)]'''
    
    
    ran1 = [k for k in range(0, n+1)] # = [0,1,....n]
    ran2 = [k for k in range(0, m+1)] # = [0,1,....m]
    
    combo = []
    
    for i in ran1:
        for j in ran2:
            combo.append((i,j))
    
    return combo
    
def generate_subset(n, num):
    '''given range [n] = [0,1,2,...n], returns the list of subsets
    of size n of [n]'''
    return list(it.combinations(range(n+1), num))

def generate_bijections(n):
    ''' given n, gives a list of possible 'subsets'of n, including 0 and n
    , EXCLUDING empty set. Used for generating idempotents'''
    
    lst = []
    for size in range(1,n+2):
        lst += generate_subset(n, size)
    
    return lst

def generate_bijections_2(n,m):
    ''' Given n and m, it gives bijections from [n] and [m] with [n] and [m] 
    that may have different size:
        [n] = [0, 1, 2, 3....n]
        [m] = [0, 1, 2, 3,...m]
        '''
    
    app = []
    if n > m:
        bij = generate_bijections_same(m,m)
        diff = n - m 
        for b in bij:
            arr = []
            for pair in b:
                arr.append((pair[0] + diff, pair[1]))
            app.append(arr)
    
    if n < m:
        bij = generate_bijections_same(n,n)
        diff = m - n
        for b in bij:
            arr = []
            for pair in b:
                arr.append((pair[0], pair[1] + diff))
            app.append(arr)
    
    print("++++++")
    print(bij)
    set1 = set(bij)        
    return list(set1.union(set(app)))

def generate_bijections_3(a1, b, a2):
    ''' '''
    pass
    

def generate_bijections_same(n):
    ''' Given n and m, it gives bijections from [n] and [m] with [n] and [m] size the same'''

    m = n
    ran1 = [k for k in range(0, n+1)] # = [0,1,....n]
    ran2 = [k for k in range(0, m+1)] # = [0,1,....m]
    
    bij = []
     
    for i in range(1, n+1): # i goes from 1 to n
        for set1 in it.combinations(ran1, i):
            for set2 in it.combinations(ran2, i):
                for p in it.permutations(set2):
                    bij.append(list(zip(set1,p)))
    
    return bij

def intersections( dict1, dict2, itself = False):
    ''' given dict1, and dict2, returns dictionary `intersection`that contains
     unique intersection points. 
    If itself = True, it means dict1 = dict2, and need to use a diff method
    If itself = False, we assume dict 1 and dict2 are distinct(no shared
    pairs).
    
    For example, if the intersection between (1,1):(1.5,2)(in dict1
                 and (1,2):(1.5,1) in dict2
    is (1.25,1.5) then, `intersection` will contain 
    (((1,1),(1.5,2)),((1,2),(1.5,1))) : (1.25, 1.5)
    
    i.e : The Key is the pair that causes intersection points
          The Value is the intersection points
    ''' 
    
    
    intersections = {}
    
    if itself == False:
        for key_1, value_1 in dict1.items():
            pair_1 = (key_1, value_1)
            for key_2, value_2 in dict2.items():
                pair_2 = (key_2, value_2)
                intersect = doescross(pair_1, pair_2)
                if intersect != None:
                    intersections.update({(pair_1, pair_2): intersect})
    else: # dict1 = dict2
        intersections = same_set_intersections(dict1)
    
    return intersections
                
def same_set_intersections(dict1):
    '''Given pairs by dict1, it gives the list of intersections by pairs in
    dict1 themselves.'''
    
    intersections = {}
    
    n = len(dict1)
    pairs = list(dict1.items())
    combinations = generate_subset(n-1,2) # generate combinations of [0,1,...n-1] with size 2
    
    
    for combo in combinations:
        if combo[0] == combo[1]:
            print("something is wrong, tuple has duplicate values")
        i = combo[0]
        j = combo[1]
        intersect = doescross(pairs[i],pairs[j])
        if intersect != None:
            intersections.update({(pairs[i],pairs[j]): intersect})
    
    return intersections
                
def simple_intersections( dict1, dict2, itself = False): # cb and fix positional arguments
    ''' given dict1, and dict2, returns number of unique intersection points 
    
    ''' 
    return len(intersections(dict1, dict2, itself))

def get_domain(tup):
    '''given a tuple of pairs(injective), it returns list of domains
    For example, if tup = ((1,2), (2,3), (4,2)), 
    then this method returns [1,2,4]'''
    
    dom = []
    for pair in tup:
        dom.append(pair[0])
    return dom
    
def get_range(tup):
    '''given a tuple of pairs(injective), it returns list of domains
    For example, if tup = ((1,2), (2,3), (4,4)), 
    then this method returns [2,3,4]'''
    
    ran = []
    for pair in tup:
        ran.append(pair[1])
    return ran

def get_domain_dict(dic):
    ''' either takes a dictionary object or a list of dictionary objects
    and return the domain.'''
    print("+++++ENTERING DOMAIN ++++++")
    dom = set()
    if isinstance(dic,list):
        for cur_dic in dic:
            dom.update(set(cur_dic.keys()))
    else:
        dom = set(dic.keys())

    print(dom)
    return dom

def get_range_dict(dic):
    ''' either takes a dictionary object or a list of dictionary objects
    and return the range.'''
    
    print("+++++++ ENTERING RANGE++++++++")
    ran= set()
    if isinstance(dic,list):
        for cur_dic in dic:
#            print(set(cur_dic.values()))
#            print("current ran =")
#            print(ran)
#            print("++++")
            ran.update(set(cur_dic.values()))
    else:
        ran = set(dic.values())
    print(ran)
    return ran
def is_bijection(tup):
    ''' given a tuple of pairs, it returns whether the given tuple is 
    is an idempotent
    For example, 
    1. if tup = ((1,2),(2,1)) then returns false
    2. if tup = ((1,1), (2,2)) then returns true
    '''

def safeMultiply(a, b):
    """Safely multiply the two sides using __mul__ and __rmul__. Return
    NotImplemented if both fails.
    """
    try:
        prod = a.__mul__(b)
    except TypeError:
        prod = NotImplemented
    if prod is NotImplemented:
        try:
            prod = b.__rmul__(a)
        except TypeError:
            prod = NotImplemented
    return prod

# Constants for left and right action
ACTION_LEFT, ACTION_RIGHT = 0, 1
def sideStr(side):
    if side == ACTION_LEFT: return "LEFT"
    else: return "RIGHT"
    
#
#print(orientation_i(((1,3),(2,2)),False))
#print("----1")
#print(orientation_i(((1,2),(2,3)),False))
#print("----2")
#print(orientation_i(((0,3),(1,3)),True))
#print("----3")
#print(orientation_i(((0,2),(1,2)),True))
#print("-----4")
#print(orientation_i(((2,1),(2,0)),False))
#print("----5")
#
#

#print("--------------")
#print(orientation(((0,0),(1,1)),True))
#print("--------------")
#print(orientation(((0,0),(0,1)),True))
#print("----")
#print(orientation(((0,0),(0,1)),False))
#print("----")
#print(orientation_i(((0,0),(0,1)),True))
#print("----")
#print(orientation_i(((0,0),(0,1)),False))

#print("----dictionary test starts.----")
#dict1=    {
#  "x": 0,
#  "y": 1,
#  "z": 2
#}
#
#dict2=    {
#  "x": 1,
#  "y": 1,
#  "z": 2
#}


#test code
#
##print(__add__(dict1, [dict2, dict2]))  #has error
#print(__iadd__(dict1,dict2))
##print(__sub__(dict1,[dict2,dict2]))
#print(__isub__(dict1,dict2))
#print(accumulate(dict1,[dict2, dict2]))
#print(_dictMult(dict2, dict2))
#
#print(iter(dict1).__next__())
##F2 = ModNRing(2)
##ZZ = Integer()

#print({1,3,4} not in {1,2,3,4,5})
#a = (1,2,3,4,5)
#b = (1,3,4)
#print(complement(a,b))

dict1= SummableDict({
  "x": 1,
  "y": 2,
  "z": 1
})

dict2=   SummableDict({
  "x": 1,
  "y": 4,
  "z": 1
})

dict3=   SummableDict({
  "x": 1,
  "y": 4,
  "z": 1
})

#a = generate_bijections_2(2,3)
#print(a)
##
#print(type([dict2, dict2]))
#
## Note that add and sub is only for two, and i means 'other' is a list 
#print("___add___")
#print(__add__(dict1,dict2))
#print("___iadd___")
#print(__iadd__(dict1,[dict2, dict3]))
#print("___sub___")
#print(__sub__(dict1,dict2))
#print("____isub____")
#print(__isub__(dict1,dict2))
#
##a = generate_subset(5,2)
##print(a)
##
#dict11 = {(1,1):(1.5,2),(1,2):(1.5,1), (1,3):(1.5,3)}
#i = intersections(dict11, dict11, True)
#print(i)
#
#print(doescross_simple((4, 3), (3, 4)))

#lst = generate_bijections(3)
#
#print(lst)
#print(len(lst))

#print(generate_bijections(2))
#print(generate_bijections_2(1,2))
#print(len(generate_bijections_2(1,2)))

print(len(combinations(2, 3)))
