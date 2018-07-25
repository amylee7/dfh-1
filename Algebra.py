#from numbers import Number
#import heapq
#from numbers import Number
from utility import NamedObject, SummableDict,safeMultiply
#from utility import fracToInt, memorize, memorizeHash, safeMultiply
#from utility import F2
from numbers import Number

class FreeModule:
    '''Represents a free module over some ring.'''
    
    def __init__(self,ring):
        '''Specifies the ring. Ring should be a subclass of Ring or a python number'''
        self.ring = ring
    
    def getGenerators(self):
        '''returns a list of all generators. Need not be implemented for every module.'''
        
        raise NotImplementedError("Method getGenerators not implemented.")

class Generator:
    '''Represents a generator of a free module.
    Implement __eq__, __ne__, and __hash__.'''
    
    def __init__(self,parent):
        '''every Generator needs the info of its "parent" module'''
        self.parent = parent
    
    def elt(self, coeff = 1):
        '''Returns the element coeff * self'''
        return self.ELT_CLASS({self: coeff})  # returns the element(summable dict)
                                        # object with key self, value coeff
    
    def diff(self):
        '''returns the differential of this generator.
        It calls the diff() method of parent module.'''
        
        return self.parent.diff(self)
    
    def factor(self):
        ''' Finding all ways to factor this generator into a product of
        two generators. Calls factor() method of its parent module'''
        
        return self.parent.factor(self)

    def delta(self):
        '''Returns the delta of this generator. 
        It calls the delta() method of parent module.'''
        
        return self.parent.delta(self)
    
    def __mult__(self, other):
        '''multiples this generator by ``other`` on the left.
        That is, the right module action'. Usually returns an element
        rather than a generator''' # left?
        
        if isinstance(other,Number):
            return self.elt(other)
        elif hasattr(self.parent, "multiply"):
            return self.parent.multiply(self,other)
        else:
            return NotImplemented
    
    def __rmul__(self,other):
        '''multiplies this generator by ``other`` on the right. 
        That is, the left module ation
        Usually represents a left module action''' # right?
        
        if isinstance(other, Number):
            return self.elt(other)
        elif hasattr(self.parent, "rmultiply"):
            return self.parent.rmultiply(self, other)
        else:
            return NotImplemented


    
    def toSimpleGenerator(self, name): # ? Necessary? 
        """Convert to a SimpleGenerator with the given name. All fields are
        preserved, except ``name`` which is overwritten, and _hash_val which is
        removed, if present.
        """
        new_obj = SimpleGenerator(self.parent, name)
        new_obj.__dict__.update(self.__dict__)
        new_obj.name = name # to make sure original name is overwritten
        if hasattr(new_obj, '_hash_val'):
            del new_obj._hash_val # reset hash value
        return new_obj

class SimpleGenerator(Generator, NamedObject):
    '''Each Generator has a name. Distinguished by name.'''
    def __init__(self, parent, name):
        """Specifies name and parent module."""
        Generator.__init__(self, parent)
        NamedObject.__init__(self, name) # ?
    
class Element(SummableDict):
    '''Represents an element of a free module, as a dictionary from generators
    to coefficients. ex) a+3b will be written as {a:1,b:3}'''
    
    def __init__(self,data = None):
    
        if data is None:
            data = {}
        
        SummableDict.__init__(self,data)
        
        if self: # if the object is instantialized
            convert = self.getElt().parent.ring.convert #?
            for key,value in self.items():
                self[key] = convert(value)
    
    def __str__(self):
        #string method for Element Class
        if self == 0:
            return "0"
        terms = []
        for gen, coeff in self.items():
            if coeff == 1:
                terms.append(str(gen))
            else:
                terms.append(str(coeff)+"*"+str(gen))
        return "+".join(terms)
    
    def __repr__(self):
        return str(self)
    
    def __mul__(self, other):
        # First try multiplying each coefficient with other, using the function
        # in SummableDict.
        result = SummableDict.__mul__(self, other)
        if result != NotImplemented:
            return result

        # Now try to multiply each key by other on the right.
        result = E0
        for k, v in self.items():
            prod = safeMultiply(k, other)
            if prod is NotImplemented:
                return NotImplemented
            result += [term * (v * coeff) for term, coeff in prod.items()]
        return result
    
    
    def __rmul__(self, other):
        # First try multiplying each coefficient with other, using the function
        # in SummableDict.
        result = SummableDict.__rmul__(self, other)
        if result != NotImplemented:
            return result

        # Now try to multiply key by other on the left.
        result = E0
        for k, v in self.items():
            prod = safeMultiply(other, k)
            if prod is NotImplemented:
                return NotImplemented
            result += [term * (v * coeff) for term, coeff in prod.items()]
        return result
    
    def diff(self):
        ''' Returns the differential of this element.'''
        return sum([coeff * gen.diff() for gen, coeff in self.items()], E0)
  
# Name of the class for elements containing this generator
Generator.ELT_CLASS = Element
# Short-hand for empty element
E0 = Element()

class ChainComplex(FreeModule): # ask akram : why is this extension of free module
    '''Represents a general chain complex. '''

    def diff(self, gen):
        ''' Returns the differential of a generator.'''
        raise NotImplementedError("Differential not implemented.")

class SimpleChainComplex(ChainComplex): 
    '''Represents a chain complex with a finite number of generators, with 
    explicitly stored generating set and differential.Generating set is stored
    as a python set, and differential is stored as a dictionary mapping
    from generators to elements. Each generator must be a key in the dictionary
    (even if its differential is zero)'''
    
    def __init__(self,ring):
        ''' Initialize an empty chain complex.'''
        ChainComplex.__init__(self,ring)
        self.generators = set()
        self.differential = dict()
    
    def __str__(self):
        result = "Chain complex. \n"
        for k,v in self.differential.items():
            result += "d(%s) = %s \n" %(k,v)
        return result

    def __repr__(self):
        return str(self)
    
    def __len__(self):
        return len(self.generators)
    
    def diff(self, generator):
        return self.differential[generator]
    
    def diffElt(self, elt):
        '''Return the differential of Element elt of this SimpleChainComplex'''
        # applying linearity to element
        
        answer = E0
        for x in elt.keys():
            answer += elt[x]*self.diff(x)
        return answer
    
    def getGenerators(self):
        return list(self.generators)
    
    def reindex(self): # cb and fill in
        '''Replace the generators by simple generators indexed by integers.
        The names of the new generators are 'g1', 'g2', etc.'''
        
        gen_list = list(self.generators)
        new_gen_list = []
        
        # Dictionary mapping original generators to new ones
        translate_dict = dict()
        for i in range(len(gen_list)):
            new_gen = gen_list[i].toSimpleGenerator("g%d"%(i+1))
            new_gen_list.append(new_gen)
            translate_dict[gen_list[i]] = new_gen
        self.generators = set(new_gen_list)
        new_diff = dict()
        for gen, dgen in self.differential.items():
            new_diff[translate_dict[gen]] = dgen.translateKey(translate_dict)
        self.differential = new_diff
        
        #translate Maslow grading
        if hasattr(self, "m_grading"):
            new_grading = dict()
            for gen, gr in self.m_grading.items():
                if gen in translate_dict: # gen is still in chain complex
                    new_grading[translate_dict[gen]] = gr
            self.m_grading = new_grading
        
        #translate Alexander grading
        if hasattr(self, "a_grading"):
            new_grading = dict()
            for gen, gr in self.a_grading.items():
                if gen in translate_dict: # gen is still in chain complex
                    new_grading[translate_dict[gen]] = gr
            self.a_grading = new_grading

    def addGenerator(self, generator):
        '''Add an generator. No effect if the generator already exists.'''
        
        assert generator.parent == self
        self.generators.add(generator)
        #if not self.differential.has_key(generator): # has key removed in python 3
        if generator not in self.differential:
            self.differential[generator] = E0    
    
    def addDifferential(self, gen_from, gen_to, coeff):
        ''' Add coeff * gen_to to the differential of gen_from.
        i.e diff(gen_from) = coeff * gen_to
        gen_from and gen_to must be generators of this complex
        '''
        
        assert gen_from.parent == self and gen_to.parent == self
        self.differential[gen_from] += coeff * gen_to
        
    #cb and recode
    def simplify(self, find_homology_basis = False, cancellation_constraint = None): 
        pass
    
    def checkDifferential(self):
        '''Check the relations d^2 for differentials.'''
        for gen in self.generators:
            assert gen.diff().diff() == 0
    
    def checkGrading(self): # cb and fix grading[x] - where is the grading method?
        ''' Check if the grading is consistent with differentials
        i.e., Alexander is preserved, and Maslov down by one '''
        for x in self.generators:
            for y, coeff in x.diff().items():
                assert self.m_grading[x] - 1 == self.m_grading[y] \
                        and self.a_grading[x] == self.a_grading[y]
    
    def getGradingInfo(self):
        ''' Shows the distribution of gradings in an easy-to-read format.'''
        
        # cb and fill in - diff
        
    def copy(self, returnDicts = False): # cb and fill in 
        pass
    
    def id(self): # cb and remove if unnecessary
        pass

#class SimpleChainMorphism: # cb and remove if uncessary
#    '''Represents a morphism between two simple chain complexes (that may not be distinct)
#    Need not be chain map . Represented explicitly.'''
    

class DGAlgebra(ChainComplex):
    '''Represents a general differential-graded algebra.'''
    def multiply(self, gen1, gen2):
        '''Returns the product of gen1 and gen 2, as an algebra element.'''
        raise NotImplementedError("multiply not implemented.")
        
class Tensor(FreeModule, tuple):
    '''Represents a free module whose generating set is the product of the generating set of two sides.''' 
    #two sides ?  ask akram
    #tuples
    
    def __init__(self,data):
        '''Specifies the left and right module'''
        for d in data[1:]:
            assert d.ring == data[0].ring # raise assertion if the base ring is not the same
    #Note that the tuple initialization is automatic 
        FreeModule.__init__(self,data[0].ring) # specifies the generator and the rings 

class TensorGenerator(Generator, tuple): #?
    '''Represents a generator of a free module that is a tensor product two more free modules. 
    Works as a tuple of the components'''

class TensorElement(Element):
    '''Represents an element of the tensor product of two or more modules.'''

TensorGenerator.ELT_CLASS = TensorElement

def expandTensor(prod, parent = None):
    '''Produces the tensor element formed by the tensor product of either
    generators or elements.
    
    ``prod`` is a tuple of either Generator or Element, corresponding to the
    components of the tensor product.
    
    For example, ((1*A+1*B),C) expands into 1*(A,C)+1*(B,C), and
    ((1*A-1*B),(1*C-1*D)) expands into 1*(A,C)-1*(B,C)-1*(A,D)+1*(B,D).
    
    ``parent`` specifies the Tensor module. If it is set to None, the default
    (with no additional operations defined) will be used (during the
    initialization of TensorElement).'''

def TensorDGAlgebra(Tensor, DGAlgebra):
    '''Tensor Product of DGAlgebra(Differential Graded Algebra) is a DG Algebra'''

def TensorIdempotent(tuple):
    '''Serves as idempotent to a tensor product of algebras.'''

class TensorStarGenerator(Generator, tuple):
    '''Represents a generator of the tensor star algebra - a tuple (possibly
    with zero components) of elements in the same algebra.'''

def simplifyComplex(arrows, default_coeff = 0, find_homology_basis = False,
                    cancellation_constraint = None):
    '''Simplify complex using the cancellation lemma.'''

def findRankOverF2(num_row, num_col, entries):
    '''Find rank of a matrix over F2 with the given number of rows and columns.
    entries is a list of pairs (i, j) with 0 <= i < num_row and 0 <= j < num_col
    specifying where the matrix has 1's.'''
