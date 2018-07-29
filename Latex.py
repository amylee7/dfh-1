''' Colletion of Latex Printing Code.'''

def beginDoc():
    return "\\documentclass{article}\n" + \
        "\\usepackage{tikz}\n" + \
        "\\begin{document}\n"

def endDoc():
    return "\\end{document}\n"

def beginTikz(scale):
    return "\\begin{tikzpicture} [x=%dpt, y=%dpt," % (3 * scale, scale) + \
        "baseline=(current bounding box.center)]\n"

def endTikz():
    return "\\end{tikzpicture}\n"


def showStrandAlgebraElement(sa_elt):
    ''' `sa_elt` is a local strandalgebra element to be displayed.'''
    result = ""
    tangle = sa_elt.tangle
    strands = sa_elt.strands
    
    # Header
    result += beginTikz(30)
    if sa_elt.is_left:
        for start,end in tangle.orient_left_rhalf.items():
            result += "\draw [->] [ultra thick] [color=orange] %s to [out=180, in=0] %s;\n" \
                    %(start,end)
        for start,end in tangle.orient_right_rhalf.items():
            result += "\draw [->] [ultra thick] [color=orange] %s to [out=0, in=180] %s;\n" \
                    %(start,end)   
        for start,end in strands.right_converted.items():
            result += "\draw  [thick] [color = black] %s to[out=0, in=180] %s; \n" %(start,end)
    else:
        for start,end in tangle.orient_left_lhalf.items():
            result += "\draw [->] [ultra thick] [color=orange] %s to [out=180, in=0] %s;\n" \
                %(start,end)
        for start,end in tangle.orient_right_lhalf.items():
            result += "\draw [->] [ultra thick]  [color=orange] %s to [out=0, in=180] %s;\n" \
                %(start,end)
        for start,end in strands.left_converted.items():
            result += "\draw [thick] [color = black] %s to [out=0, in=180] %s; \n" %(start,end)

    # Footer
    result += endTikz()
    
    return result

def showTensor(sd,sa_elt):
    '''shows tensor between strand diagram element and strand algebra element'''
    result = ""
    result += "\\begin{equation}\n"
    result += showStrandDiagram(sd)
    result += "~\\otimes~ \n"  
    result += showStrandAlgebraElement(sa_elt)
    result += "\\end{equation}\n"
    return result

def showStrandDiagram(sd):
    ''' `sd` is the local strand diagram to be displayed. '''
    
    result = ""
    tangle = sd.tangle
    strands = sd.strands
    
    # Header
    result += beginTikz(30) 
    # Display tangles for Strand Diagram
    for start,end in tangle.orient_left_lhalf.items():
        result += "\draw [->] [ultra thick] [color=orange] %s to [out=180, in=0] %s;\n" \
                %(start,end)
    for start,end in tangle.orient_right_lhalf.items():
        result += "\draw [->] [ultra thick]  [color=orange] %s to [out=0, in=180] %s;\n" \
                %(start,end)
    for start,end in tangle.orient_left_rhalf.items():
        result += "\draw [->] [ultra thick] [color=orange] %s to [out=180, in=0] %s;\n" \
                %(start,end)
    for start,end in tangle.orient_right_rhalf.items():
        result += "\draw [->] [ultra thick] [color=orange] %s to [out=0, in=180] %s;\n" \
                %(start,end)           
    #Display Strands----------
    for start,end in strands.left_converted.items():
        result += "\draw [thick] [color = black] %s to [out=0, in=180] %s; \n" %(start,end)
    for start,end in strands.right_converted.items():
        result += "\draw  [thick] [color = black] %s to[out=0, in=180] %s; \n" %(start,end)
        
    # Footer
    result += endTikz()
    return result
#    
#def showDAGenerator():
#    '''Show a generator of the type DA Structure.'''
#    pass
#
#def showDAStructure(dastr):
#    pass
    
