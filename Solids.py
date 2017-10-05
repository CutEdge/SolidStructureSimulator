"""
Solids - A program to analyze complex 3D systems of solid links and nodes.
==========================================================================

Update from previous drafts:
-Removed the example program.
-Placed a simpler example program into the code, for troubleshooting.
-This example program succeeds. Here's how I fixed it:
        *Remove all "f0" and "m0" connections.
        *strip the coefficients matrix K of all rows/columns that are pure zeros.
        *define the result matrix R in terms of the number of rows, not the number of columns.

"""
#Cautions:
"""============================
-No link may twist more than 10 degrees in any direction.

"""      
#Notes:
"""============================
-This version contains fewer notes. See previous versions for elaboration on any unclear details.
-"Gradient" does not mean the vector calculus term; it's the slope of a torqued or bent link.
-all vectors are defined as positive up, right, and out of the page.
-"force", "moment", "gradient", and "displacement" are all properties of each end of a link; they are not properties of a node.
        The exception is for Applied Constraints.
-Applied Constraints are termed "fn0", "mn0", "tn0", "dn0", "fn1", "mn1" ... etc.
        So, at node0, with nodal force "fn0", the static force equation might be:
                -Fl1nb - Fl2nb - Fl3na + fn0 = 0
        where -Fl1nb means "reaction force from node b of link 1"
-Each Applied constraint will have its own *separate* equation in which the "fn0" term is defined. The coefficient matrix
        will be an identity matrix, and the right side will have some vector; maybe like this: I*fn0 = [1;3;2]
-a Connector governs the interface between a node and a link.
-The 4 types of matrix are "Link", "Node", "Connector", and "Applied Constraint".

"""
#Assumptions:
"""============================
* Forces on a link sum to zero.
* Forces on a node sum to zero.
* Force flips if transferred from a link to a node: F_link = -F_node.
* Connectors of type "f0" cause forces not to transfer from a link to a node.
* Applied Force variable * 1 = Applied Force value(definition)

* Moments on a link sum to zero.
* Moments connected to a node sum to zero.
* Moment flips if transferred from a link to a node: M_link = -M_node.
* Connectors of type "m0" cause moments not to transfer from a link to a node.
* Applied Moment variable * 1 = Applied Moment value (definition)

* Gradients on a link: t = integral(M,dr)/(EI) for a given link.
* Gradients connected to a given node are all equal.
* Gradient does not flip if transferred from a link to a node: t_link = t_node.
* Connectors of type "tx" cause gradients to transfer from a link to a node.
* Applied Gradient variable * 1 = Applied Gradient value (definition)

* Displacements on a link: d = integral(t,dr) for a given link.
* Displacements connected to a given node are all equal.
* Displacement does not flip if transferred from a link to a node: d_link = d_node.
* Connectors of type "dx" cause displacements to transfer from a link to a node.
* Applied Displacement variable * 1 = Applied Displacement value (definition)

* Applied loads must consist of (Force or Displacement) and (Moment or Gradient) at every node.

"""
#Equations for a given link:
"""============================
> Fa + Fb = 0
> Ma + Mb + (r x Fb) = 0
> tb = -|r|*(r^ x Ma x r^)/(EI)   -   |r|*(r x Fb)/(2EI)   -   |r|*(r^ . Ma . r^)/(GJ)   +   ta
> db = -|r|*(r^ x Ma x r^) x r/(2EI)    -    |r|*(r x Fb x r)/(6EI)    +    (ta x r)    +    (r^ . Fb . r)/(E*A)    +    da
------->note: (r^ x Ma x r^) x r can be expressed in matrix form as cross(r,transverse(r)). For some reason, crossing matrices flips the sign.
"""
#Equations for a given node:
"""============================
> F0 + F1 + F2 + ... + F_from_each_link - f_applied_at_node = 0
> M0 + M1 + M2 + ... + M_from_each_link - m_applied_at_node = 0
"""
#Equations for a connector at a given end "a" of a given link "l", to a given node "n":
"""============================
> if type "f0": F_la = 0
> if type "m0": M_la = 0
> if type "tx": t_la - t_n = 0
> if type "dx": d_la - d_n = 0
"""
#Equations for an applied constraint "A" on a given node "n":
"""============================
> F_n = A_f
> M_n = A_m
> t_n = A_t
> d_n = A_d
"""
import numpy
from numpy import *
from numpy.linalg import norm, solve
from math import log10

def A(D): #Area
    if len(D)==1: #Assume a diameter is given
        A = pi*D[0]**2/4.
    elif len(D)==2: #Assume [base,height] of a rectangle is given
        A = D[0]*D[1]
    return A

def I(D): #Area Moment of Inertia
    if len(D)==1: #Assume a diameter is given
        I = pi*D[0]**4/64.
    elif len(D)==2: #Assume [base,height] of a rectangle is given
        I = 1/12.*D[0]*D[1]**3
    return I

def J(D): #Polar Moment of Inertia
    if len(D)==1: #Assume a diameter is given
        J = 2.*pi*D[0]**4/64.
    elif len(D)==2: #Assume [base,height] of a rectangle is given
        J = 1/12.* ( D[0]*D[1]**3 + D[1]*D[0]**3 )
    return J

def unit(v):
        #divides a vector or matrix by its own magnitude.
        return v/norm(v)

def parallel(vector,direction):
        #projection of 'vector' parallel to 'direction'
        #only works if both 'vector' and 'direction' are fully defined.
        return float(dot(vector,direction))/dot(direction,direction)*direction

def perpendicular(vector,direction):
        #projection of 'vector' perpendicular to 'direction' in the common plane.
        #only works if both 'vector' and 'direction' are fully defined.
        return vector - parallel(vector,direction)

def rotate(vector,axis,degrees,digits=10):
        vector = array(vector)
        axis = array(axis)
        angle = radians(degrees)
        if norm(vector)==0:
                rotatedVector = zeros([1,3])
        elif (vector==axis).all():
                rotatedVector = array([vector])
        else:
                vpp = perpendicular(vector,axis)
                vpl = parallel(vector,axis)
                w = cross(axis,vpp)
                rotatedVector = vpl + round(cos(angle),digits)*vpp + round(sin(angle),digits)*unit(w)*norm(vpp)
        return rotatedVector
        

def matrixCross(v):
        #converts a 3-vector into a matrix. This matrix can be multiplied
        #by a 3x1 matrix to cross the vector with that matrix.
        #[v x ?] = matrixCross(v) * ?
	try: m = matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
	except:
		v = v[0]
		m = matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
	return m

def doubleCross(v):
        #[v x ? x v] = doubleCross(v) * ?
        #gives the projection of ? perpendicular to ? in the shared plane, times v^2.
        #in matrix form.
        return ( array(matrix(v) * matrix(v).T)[0][0] * matrix(eye(3)) - matrix(v).T * matrix(v) )

def matrixDot(v):
        #[v . v] = matrixDot(v)
        return matrix(v) * matrix(v).T

def doubleDot(v):
        #[v . ? . v] = doubleDot(v) * ?
        #gives the projection of ? parallel to v, times v^2.
        #in matrix form.
        return matrix(v).T * matrix(v)

def axial(v):
        #[v^ . ? . v^] = axial(v) * ?
        #gives the projection of ? parallel to v.
        #in matrix form.
        return doubleDot(v)/norm(v)**2

def transverse(v):
        #[v^ x ? x v^] = transverse(v) * ?
        #returns the projection of some unknown vector '?', perpendicular to the vector 'v',
        #in the plane shared by 'v' and '?'
        #expresses the answer in matrix form.
        return doubleCross(v)/norm(v)**2

def issame(list1,list2):
        for eachitem in (list1==list2):
                if eachitem:
                        answer = True
                else:
                        answer = False
                        break
        return answer

def saveMatrix(M,filename,delimeter = "\t"):
        f = file(filename,"w")
        for i in M:
                for j in i:
                        f.write(str(j) + delimeter)
                f.write("\n")
        f.close()

def stripZeroRows(a):
        global rowsRemoved
        global rowsKept
        rowsRemoved=[]
        h = shape(a)[0]
        w = shape(a)[1]
        rowsKept=range(h)
        i=0
        for i in range(h):
                if (a[i]==zeros([1,w])).all():
                        rowsRemoved.append(i)
                        rowsKept.remove(i)
        return a[rowsKept]

def stripZeroCols(a):
        global colsRemoved
        global colsKept
        colsRemoved=[]
        h = shape(a)[0]
        w = shape(a)[1]
        colsKept = range(w)
        j=0
        for j in range(w):
                if (a.T[j]==zeros([1,h])).all():
                       colsRemoved.append(j)
                       colsKept.remove(j)
        return a.T[colsKept].T

def stripZeroRC(a):
        return stripZeroCols(stripZeroRows(a))

def smallestRow(a):
        return matrix([norm(a[0]),norm(a[1]),norm(a[2])]).T.argmin()
        

class Connection:
        #at any connection, some value, in some direction, is set to either zero, or "the value at that node in that direction".
        #       forces and moments will be set to zero.
        #       thetas and deltas will be set to "the value at that node in that direction". (negative projected matrix positioned in the column of the node in question).
        #
        #some examples:
        #       3d Pin (U-Joint):       m = [0,0,0] d = d_node
        #       Cylindrical Guide:      f_axial = [0,0,0], m_axial = [0,0,0], t_transverse = t_node_trans, d_transverse = d_node_trans
        #       Rectangular Guide:      f_axial = [0,0,0], t_axial = t_node_axial, t_transverse =  t_node_trans, d_transverse =  d_node_trans
        #       2d Pin (Hinge):         m_pin_axis = [0,0,0], t_pin_axis = t_node_pin_axis, d = d_node
        #       Ball on Wall:           f_planar = [0,0,0], m = [0,0,0], d_normal = d_node_normal
        #
        #some terminology:
        #       f_axial         the projection of the force vector at this node's end of the link in question, along the main axis of the link.
        #       f_transverse    the projection of the force vector at this node's end of the link in question, perpendicular to the axis of the link.
        #       f_planar        the same as f_transverse, but referring to the force vector proejcted to the plane of a wall or surface.
        #       f_normal        the same as f_axial, but referring to the normal vector of a plane.
        #
        #the user must pass a full 3x3 matrix to the "projection" argument. The user must decide before-hand whether this connection is "axial" or "planar" or whatever,
        #       and use the axial() or transverse(), or some other function, to generate a projection matrix that can later be multiplied
        #       matrix-style with an unknown (3x1) matrix.
        #
        #For instance, if I wanted to create a Cylindrical Guide on link0 to node0:
        #       nodes[0].addConnection(links[0],"f",axial(links[0].r))
        #       nodes[0].addConnection(links[0],"m",axial(links[0].r))
        #       nodes[0].addConnection(links[0],"t",transverse(links[0].r))
        #       nodes[0].addConnection(links[0],"d",transverse(links[0].r))
        
        def __init__(self,node,link,xtype,project):
                self.node = node
                self.link = link
                self.xtype = xtype
                self.project = project
                self.projection = eye(3)
                if xtype == "f0": pass
                elif xtype == "m0": pass
                elif xtype == "tx": pass
                elif xtype == "dx": pass
                else: print "ERROR: Connection between node ",self.node," and link ",self.link," has a bad type. Type must be f0, m0, tx, or dx."
                if project == "all":pass
                elif project == "axial":pass
                elif project == "transverse":pass
                else: print "ERROR: Connection between node ",self.node," and link ",self.link," has a bad projection. 'project' must equal 'all', 'axial', or 'transverse'."
                return None
        def show(self):
                print self.node, self.link, self.xtype, self.project,"\n", self.projection

class Node:
        def __init__(self,i,r,constraint=None,F=[0,0,0],M=[0,0,0],t=[None,None,None],d=[None,None,None]):
                self.i = int(i) #index of this node
                if constraint == "fix":
                    t=[0,0,0]
                    d=[0,0,0]
                elif constraint == "pin":
                    t=[None,None,None]
                    d=[0,0,0]
                elif constraint == None or constraint == "None":
                    pass
                else:
                    print "ERROR: Node ",self.i," has an invalid constraint. Must be 'fix', 'pin', or None."
                    
                self.r = array(r) #global location of this node
                self.constraint = constraint
                self.F = c_[F] #any applied force on this node
                self.M = c_[M] #any applied force on this node
                self.t = c_[t] #any gradient constraint on this node
                self.d = c_[d] #any displacement constraint on this node
                self.connections = []
                return None
        def show(self):
                print "i:", self.i
                print "r:", self.r
                print "constraint:", self.constraint
                print "F:", self.F
                print "M:", self.M
                print "t:", self.t
                print "d:", self.d
                print "\n"
                return None

class Link:
        def __init__(self,i,na,nb,A,I,J,E,G):
                self.i = int(i)   #index of this link
                self.na = int(na) #index of one node
                self.nb = int(nb) #index of other node
                self.A = float(A) #Area of cross section of link
                self.I = float(I) #Area Moment of Inertia
                self.J = float(J) #Polar Area Moment of Inertia
                self.E = float(E) #Elastic Modulus
                self.G = float(G) #Shear Modulus

                self.r = zeros(3) #initialize an empty vector to store the length vector.
                return None
        def show(self):
                print "i:", self.i
                print "na:", self.na
                print "nb:", self.nb
                print "A:", self.A
                print "I:", self.I
                print "J:", self.J
                print "E:", self.E
                print "G:", self.G
                print "\n"
                return None
        
class System:
        def __init__(self,nodes,links,connections=[],verbose=False):
                self.nodes=nodes #list of nodes
                self.links=links #list of links
                self.connections=connections #list of connections. Connections must be assigned *after* the system is created, so it will have already solved the link vectors.
                self.verbose=verbose
                self.makeVectors()
                return None

        def makeVectors(self):
                for l in self.links:
                        l.r = self.nodes[l.nb].r - self.nodes[l.na].r #define r, the vector from na to nb, for each link.
                        if self.verbose: print "Link ",l.i," vector is ",l.r
                return None

        def testLinks(self):
                #Test all links:
                for l in self.links:
                        if l.na == l.nb: print "ERROR: Link ",l.i," is circular: both ends touch the same node."
                return None

        def makeConnections(self):
                if self.verbose: print "Making connections..."
                #Test all connections:
                i = 0
                for c in self.connections:
                        thislink = self.links[c.link]
                        thisnode = self.nodes[c.node]
                        if not (thislink.na == thisnode.i or thislink.nb == thisnode.i):
                                print "ERROR: Connection ",i," failed: neither end of link ",c.link," equals node ",c.node,"."
                        i+=1
                        
                        #assign projection matrices to each connection.
                        #connections created during this process will be appended to the end of the list, so they will all be assigned too.
                        a = axial(thislink.r)
                        t = transverse(thislink.r)
                        srt = smallestRow(t)
                        if c.project == "all":
                                c.projection = eye(3)
                        elif c.project == "axial":
                                c.projection = zeros([3,3])
                                c.projection[srt] = a[srt]
                                if c.xtype == "tx": #if an axial gradient is defined, create a zero transverse moment
                                        self.connections.append(Connection(thisnode.i,thislink.i,"m0","transverse"))
                                elif c.xtype == "dx": #if an axial displacement is defined, create a zero transverse force
                                        self.connections.append(Connection(thisnode.i,thislink.i,"f0","transverse"))
                        elif c.project == "transverse":                                
                                c.projection = t.copy()
                                c.projection[srt] = zeros([1,3])
                                if c.xtype == "tx": #if a transverse gradient is defined, create a zero axial moment
                                        self.connections.append(Connection(thisnode.i,thislink.i,"m0","axial"))
                                elif c.xtype == "dx": #if a transverse displacement is defined, create a zero axial force
                                        self.connections.append(Connection(thisnode.i,thislink.i,"f0","axial"))

                nc = len(self.connections)
                nl = len(self.links)
                if nc != 4*nl:
                        print "ERROR: Must have 4 connections per link. ",nc," connections exist, and ",nl," links exist."
                        print "         (unless there are partial connections, where one or two rows are zero... then more connections are needed)."
                        print "         If there is a link-end with no tx, it needs an m0; if there is a link-end with no dx, it needs an f0."
                        print "         Connections are shown below:"
                        for l in self.links:
                                print "Connections on Link ",l.i,":"
                                for c in self.connections:
                                        if c.link==l.i: c.show()
                                print "\n\n"
                else: print "   ",nc," connections and ",nl," links."
                return None

        def allocate(self):

                #Allocate space for all mini-matrices:
                self.I3 = eye(3)
                self.I3l = eye(3*len(self.links))
                self.I3n = eye(3*len(self.nodes))
                self.Znl = zeros([3*len(self.nodes),3*len(self.links)])
                self.Zln = zeros([3*len(self.links),3*len(self.nodes)])

                #Allocate space for the total combined matrix
                self.K = 0*eye(12*3*len(self.links)+4*3*len(self.nodes))
                self.R = zeros([12*3*len(self.links)+4*3*len(self.nodes),1])
                self.names = zeros([12*3*len(self.links)+4*3*len(self.nodes),1]).astype("S100")
                
                #Coefficients of forces on links, at node a and node b:
                self.Kf_lna = self.I3l.copy()
                self.Kf_lnb = self.I3l.copy()

                #Coefficients of forces on nodes:
                self.Kf_nna = self.Znl.copy() #initiate now, fill later with carefully placed I3 matrices.
                self.Kf_nnb = self.Znl.copy() #initiate now, fill later with carefully placed I3 matrices.
                self.Kfn = -self.I3n.copy()
                
                #Coefficients of moments derived from force:
                self.Kfr_lnb = 0*self.I3l.copy() #initiate now, fill later with r terms
                #Coefficients of moments on links, at node a and node b:
                self.Km_lna = self.I3l.copy()
                self.Km_lnb = self.I3l.copy()

                #Coefficients of moments on nodes:
                self.Km_nna = self.Znl.copy() #initiate now, fill later with carefully placed I3 matrices.
                self.Km_nnb = self.Znl.copy() #initiate now, fill later with carefully placed I3 matrices.
                self.Kmn = -self.I3n.copy()

                #Coefficients of thetas derived from force and moment:
                self.Kfrr_lnb = 0*self.I3l.copy() #initiate now, fill later with r terms
                self.Kmr_lna = 0*self.I3l.copy() #initiate now, fill later with r terms
                #Coefficients of thetas on links, at node a and node b:
                self.Kt_lna = -self.I3l.copy()
                self.Kt_lnb = self.I3l.copy()

                #Coefficients of deltas derived from force, moment, and gradient:
                self.Kfrrr_lnb = 0*self.I3l.copy() #initiate now, fill later with r terms
                self.Kmrr_lna = 0*self.I3l.copy() #initiate now, fill later with r terms
                self.Ktr_lna = 0*self.I3l.copy() #initiate now, fill later with r terms
                #Coefficients of deltas on links, at node a and node b:
                self.Kd_lna = -self.I3l.copy()
                self.Kd_lnb = self.I3l.copy()

                #Coefficients of Connectors:
                self.Kf_clna = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Kf_cna = self.Zln.copy() #initiate now, fill later with carefully placed projections
                self.Kf_clnb = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Kf_cnb = self.Zln.copy() #initiate now, fill later with carefully placed projections
                self.Km_clna = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Km_cna = self.Zln.copy() #initiate now, fill later with carefully placed projections
                self.Km_clnb = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Km_cnb = self.Zln.copy() #initiate now, fill later with carefully placed projections
                self.Kt_clna = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Kt_cna = self.Zln.copy() #initiate now, fill later with carefully placed projections
                self.Kt_clnb = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Kt_cnb = self.Zln.copy() #initiate now, fill later with carefully placed projections
                self.Kd_clna = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Kd_cna = self.Zln.copy() #initiate now, fill later with carefully placed projections
                self.Kd_clnb = 0*self.I3l.copy() #initiate now, fill later with projections
                self.Kd_cnb = self.Zln.copy() #initiate now, fill later with carefully placed projections

                #Coefficients of Applied Constraints:
                self.Kf_an = 0*self.I3n.copy() #initiate now, fill later with carefully placed I3 matrices. if a given node has an f coefficient, it must not have a d coefficient.
                self.Km_an = 0*self.I3n.copy() #initiate now, fill later with carefully placed I3 matrices. if a given node has an m coefficient, it must not have a t coefficient.
                self.Kt_an = 0*self.I3n.copy() #initiate now, fill later with carefully placed I3 matrices. if a given node has a t coefficient, it must not have an m coefficient.
                self.Kd_an = 0*self.I3n.copy() #initiate now, fill later with carefully placed I3 matrices. if a given node has a d coefficient, it must not have an f coefficient.

                #Applied Constraints:
                self.FD_applied = zeros([3*len(self.nodes),1]) #initiate now, fill later
                self.MT_applied = zeros([3*len(self.nodes),1]) #initiate now, fill later
                return None

        def place(self,smallMatrix,bigMatrix,location=(0,0)):
                #places a small matrix into a big matrix at a location,
                #       overwriting any values in the big matrix that are in the way.
                startRow = location[0]
                startCol = location[1]
                endRow = startRow + smallMatrix.shape[0]
                endCol = startCol + smallMatrix.shape[1]
                bigMatrix[startRow:endRow,startCol:endCol] = smallMatrix.copy()
                return bigMatrix
        
        def buildMatrix(self):
                self.testLinks()
                self.makeConnections()
                self.allocate()

                
                for L in self.links:
                        if self.verbose: print "Link ",L.i," equations..."
                        #Equations for a given link "L" (excluding pre-defined parts):
                        """
                        > Fa + Fb = 0
                        > Ma + Mb + (r x Fb) = 0
                        > tb = -|r|*(r^ x Ma x r^)/(EI)   -   |r|*(r x Fb)/(2EI)   -   |r|*(r^ . Ma . r^)/(GJ)   +   ta
                        > db = -|r|*(r^ x Ma x r^) x r/(2EI)    -    |r|*(r x Fb x r)/(6EI)    +    (ta x r)    +    |r|*(r^ . Fb . r^)/(E*A)    +    da
                        ------->note: (r^ x Ma x r^) x r can be expressed in matrix form as cross(r,transverse(r)). For some reason, crossing matrices flips the order of terms.
                        """
                        self.Kfr_lnb[ 3*L.i:3*L.i+3 , 3*L.i:3*L.i+3]  = matrixCross(L.r)                        
                        
                        self.Kfrr_lnb[ 3*L.i:3*L.i+3 , 3*L.i:3*L.i+3] = matrixCross(L.r)*norm(L.r)/(2*L.E*L.I)
                        self.Kmr_lna[ 3*L.i:3*L.i+3 , 3*L.i:3*L.i+3]  = transverse(L.r)*norm(L.r)/(L.E*L.I) + doubleDot(L.r)/(L.G*L.J)

                        self.Kfrrr_lnb[ 3*L.i:3*L.i+3 , 3*L.i:3*L.i+3]= doubleCross(L.r)*norm(L.r)/(6*L.E*L.I) - axial(L.r)*norm(L.r)/(L.E*L.A)
                        self.Kmrr_lna[ 3*L.i:3*L.i+3 , 3*L.i:3*L.i+3] = cross(L.r,transverse(L.r))*norm(L.r)/(2*L.E*L.I)
                        self.Ktr_lna[ 3*L.i:3*L.i+3 , 3*L.i:3*L.i+3]  = matrixCross(L.r)

                        
                        #Equations for a given node "N" (excluding pre-defined parts):
                        for N in self.nodes:
                                if L.na == N.i:
                                        if self.verbose: print "        Node ",N.i," found at link ",L.i," side a."
                                        self.place(self.I3,self.Kf_nna,(3*N.i,3*L.i))
                                        self.place(self.I3,self.Km_nna,(3*N.i,3*L.i))
                                        break
                        for N in self.nodes:
                                if L.nb == N.i:
                                        if self.verbose: print "        Node ",N.i," found at link ",L.i," side b."
                                        self.place(self.I3,self.Kf_nnb,(3*N.i,3*L.i))
                                        self.place(self.I3,self.Km_nnb,(3*N.i,3*L.i))
                                        break




                #Equations for a connection between a given link "L" and a given node "N".
                for C in self.connections:
                        L = self.links[C.link]
                        N = self.nodes[C.node]
                        if self.verbose: print C.xtype," connection found between link ",L.i," and node ",N.i
                        ln = 3*L.i
                        nn = 3*N.i
                        if C.xtype == "f0": #define any linear slip connections
                                if L.na == N.i:
                                        self.place(C.projection,self.Kf_clna,(ln,ln))
                                elif L.nb == N.i:
                                        self.place(C.projection,self.Kf_clnb,(ln,ln))
                        if C.xtype == "m0": #define any rotational slip connections
                                if L.na == N.i:
                                        self.place(C.projection,self.Km_clna,(ln,ln))
                                elif L.nb == N.i:
                                        self.place(C.projection,self.Km_clnb,(ln,ln))
                        if C.xtype == "tx": #define any rotational binding connections
                                if L.na == N.i:
                                        self.place(C.projection,self.Kt_clna,(ln,ln))
                                        self.place(-C.projection,self.Kt_cna,(ln,nn))
                                elif L.nb == N.i:
                                        self.place(C.projection,self.Kt_clnb,(ln,ln))
                                        self.place(-C.projection,self.Kt_cnb,(ln,nn))
                        if C.xtype == "dx": #define any linear binding connections
                                if L.na == N.i:
                                        self.place(C.projection,self.Kd_clna,(ln,ln))
                                        self.place(-C.projection,self.Kd_cna,(ln,nn))
                                elif L.nb == N.i:
                                        self.place(C.projection,self.Kd_clnb,(ln,ln))
                                        self.place(-C.projection,self.Kd_cnb,(ln,nn))

                        
                #Applied constraints on a given node "N":
                for N in self.nodes:
                        if self.verbose: print "Node ",N.i," applied constraints..."
                        for xyz in range(3):
                                #if this node has an undefined displacement, it must have a defined force:
                                if N.d[xyz][0] == None:
                                        self.FD_applied[ 3*N.i + xyz ] = N.F[xyz][0]
                                        self.Kf_an[  3*N.i + xyz ,  3*N.i + xyz ] = 1
                                #if this node has a defined displacement, it must not have a defined force:
                                else:
                                        self.FD_applied[ 3*N.i + xyz ] = N.d[xyz][0]
                                        self.Kd_an[  3*N.i + xyz ,  3*N.i + xyz ] = 1

                        for xyz in range(3):
                                #if this node has an undefined gradient, it must have a defined moment:
                                if N.t[xyz][0] == None:
                                        self.MT_applied[ 3*N.i + xyz ] = N.M[xyz][0]
                                        self.Km_an[  3*N.i + xyz ,  3*N.i + xyz ] = 1
                                #if this node has a defined gradient, it must not have a defined moment:
                                else:
                                        self.MT_applied[ 3*N.i + xyz ] = N.t[xyz][0]
                                        self.Kt_an[  3*N.i + xyz ,  3*N.i + xyz ] = 1




                #Combine all sub-matrices:
                """
                Labels to the left/top tell the number of rows/columns.
                For instance "Links" at the left means "as many rows as the number of links * 3"
                
                        Links   Links     Nodes    Links   Links     Nodes    Links   Links     Nodes    Links   Links     Nodes
                Links [Kf_lna  Kf_lnb    0        0        0        0        0        0        0        0        0        0     ]        [f_lna]   [0]
                Nodes [Kf_nna  Kf_nnb    Kfn      0        0        0        0        0        0        0        0        0     ]        [f_lnb]   [0]
                Links [0       Kfr_lnb   0        Km_lna   Km_lnb   0        0        0        0        0        0        0     ] *      [fn   ] = [0]
                Nodes [0       0         0        Km_nna   Km_nnb   Kmn      0        0        0        0        0        0     ]        [m_lna]   [0]
                Links [0       Kfrr_lnb  0        Kmr_lna  0        0        Kt_lna   Kt_lnb   0        0        0        0     ]        [m_lnb]   [0]
                Links [0       Kfrrr_lnb 0        Kmrr_lna 0        0        Ktr_lna  0        0        Kd_lna   Kd_lnb   0     ]        [mn   ]   [0]
                Nodes [0       0         Kf_an    0        0        0        0        0        0        0        0        Kd_an ]        [t_lna]   [FD_applied]
                Nodes [0       0         0        0        0        Km_an    0        0        Kt_an    0        0        0     ]        [t_lnb]   [MT_applied]
                Links [Kf_clna 0         Kf_cna   0        0        0        0        0        0        0        0        0     ]        [tn   ]   [0]
                Links [0       Kf_clnb   Kf_cnb   0        0        0        0        0        0        0        0        0     ]        [d_lna]   [0]
                Links [0       0         0        Km_clna  0        Km_cna   0        0        0        0        0        0     ]        [d_lnb]   [0]
                Links [0       0         0        0        Km_clnb  Km_cnb   0        0        0        0        0        0     ]        [dn   ]   [0]
                Links [0       0         0        0        0        0        Kt_clna  0        Kt_cna   0        0        0     ]                  [0]
                Links [0       0         0        0        0        0        0        Kt_clnb  Kt_cnb   0        0        0     ]                  [0]
                Links [0       0         0        0        0        0        0        0        0        Kd_clna  0        Kd_cna]                  [0]
                Links [0       0         0        0        0        0        0        0        0        0        Kd_clnb  Kd_cnb]                  [0]
                """                
                NL = len(self.links)*3
                NN = len(self.nodes)*3

                r = 0
                for l in self.links:
                        self.names[r + l.i*3]="Force on Link "+str(l.i)+" at node-end a."
                        self.names[r + NL + l.i*3]="Force on Link "+str(l.i)+" at node-end b."
                for n in self.nodes:
                        self.names[r + NL + NL + n.i*3]="Force on Node "+str(n.i)+"."
                r = r + NL*2 + NN
                for l in self.links:
                        self.names[r + l.i*3]="Moment on Link "+str(l.i)+" at node-end a."
                        self.names[r + NL + l.i*3]="Moment on Link "+str(l.i)+" at node-end b."
                for n in self.nodes:
                        self.names[r + NL + NL + n.i*3]="Moment on Node "+str(n.i)+"."
                r = r + NL*2 + NN
                for l in self.links:
                        self.names[r + l.i*3]="Gradient of Link "+str(l.i)+" at node-end a."
                        self.names[r + NL + l.i*3]="Gradient of Link "+str(l.i)+" at node-end b."
                for n in self.nodes:
                        self.names[r + NL + NL + n.i*3]="Gradient of Node "+str(n.i)+"."
                r = r + NL*2 + NN
                for l in self.links:
                        self.names[r + l.i*3]="Displacement of Link "+str(l.i)+" at node-end a."
                        self.names[r + NL + l.i*3]="Displacement of Link "+str(l.i)+" at node-end b."
                for n in self.nodes:
                        self.names[r + NL + NL + n.i*3]="Displacement of Node "+str(n.i)+"."
                        
                
                r = 0
                c = 0
                self.place(self.Kf_lna,self.K,(r,c))
                r = NL
                self.place(self.Kf_nna,self.K,(r,c))
                r = NL*4 + NN*4
                self.place(self.Kf_clna,self.K,(r,c))
                
                c = c + NL
                r = 0
                self.place(self.Kf_lnb,self.K,(r,c))
                r = NL
                self.place(self.Kf_nnb,self.K,(r,c))
                r = NL + NN
                self.place(self.Kfr_lnb,self.K,(r,c))
                r = 2*NL + 2*NN
                self.place(self.Kfrr_lnb,self.K,(r,c))
                r = 3*NL + 2*NN
                self.place(self.Kfrrr_lnb,self.K,(r,c))
                r = 5*NL + 4*NN
                self.place(self.Kf_clnb,self.K,(r,c))

                c = c + NL
                r = NL
                self.place(self.Kfn,self.K,(r,c))
                r = 4*NL + 2*NN
                self.place(self.Kf_an,self.K,(r,c))
                r = 4*NL + 4*NN
                self.place(self.Kf_cna,self.K,(r,c))
                r = 5*NL + 4*NN
                self.place(self.Kf_cnb,self.K,(r,c))

                c = c + NN
                r = NL + NN
                self.place(self.Km_lna,self.K,(r,c))
                r = 2*NL + NN
                self.place(self.Km_nna,self.K,(r,c))
                r = 2*NL + 2*NN
                self.place(self.Kmr_lna,self.K,(r,c))
                r = 3*NL + 2*NN
                self.place(self.Kmrr_lna,self.K,(r,c))
                r = 6*NL + 4*NN
                self.place(self.Km_clna,self.K,(r,c))

                c = c + NL
                r = NL + NN
                self.place(self.Km_lnb,self.K,(r,c))
                r = 2*NL + NN
                self.place(self.Km_nnb,self.K,(r,c))
                r = 7*NL + 4*NN
                self.place(self.Km_clnb,self.K,(r,c))

                c = c + NL
                r = 2*NL + NN
                self.place(self.Kmn,self.K,(r,c))
                r = 4*NL + 3*NN
                self.place(self.Km_an,self.K,(r,c))
                r = 6*NL + 4*NN
                self.place(self.Km_cna,self.K,(r,c))
                r = 7*NL + 4*NN
                self.place(self.Km_cnb,self.K,(r,c))

                c = c + NN
                r = 2*NL + 2*NN
                self.place(self.Kt_lna,self.K,(r,c))
                r = 3*NL + 2*NN
                self.place(self.Ktr_lna,self.K,(r,c))
                r = 8*NL + 4*NN
                self.place(self.Kt_clna,self.K,(r,c))

                c = c + NL
                r = 2*NL + 2*NN
                self.place(self.Kt_lnb,self.K,(r,c))
                r = 9*NL + 4*NN
                self.place(self.Kt_clnb,self.K,(r,c))

                c = c + NL
                r = 4*NL + 3*NN
                self.place(self.Kt_an,self.K,(r,c))
                r = 8*NL + 4*NN
                self.place(self.Kt_cna,self.K,(r,c))
                r = 9*NL + 4*NN
                self.place(self.Kt_cnb,self.K,(r,c))

                c = c + NN
                r = 3*NL + 2*NN
                self.place(self.Kd_lna,self.K,(r,c))
                r = 10*NL + 4*NN
                self.place(self.Kd_clna,self.K,(r,c))

                c = c + NL
                r = 3*NL + 2*NN
                self.place(self.Kd_lnb,self.K,(r,c))
                r = 11*NL + 4*NN
                self.place(self.Kd_clnb,self.K,(r,c))

                c = c + NL
                r = 4*NL + 2*NN
                self.place(self.Kd_an,self.K,(r,c))
                r = 10*NL + 4*NN
                self.place(self.Kd_cna,self.K,(r,c))
                r = 11*NL + 4*NN
                self.place(self.Kd_cnb,self.K,(r,c))

                c = 0
                r = 4*NL + 2*NN
                self.place(self.FD_applied,self.R,(r,c))
                r = 4*NL + 3*NN
                self.place(self.MT_applied,self.R,(r,c))
                
                self.K = stripZeroRC(self.K)
                self.names = self.names[colsKept]
                self.R = self.R[rowsKept]
                return None

        def solveMatrix(self,effective_zero = 1e-10,precision=6):
                self.solution = solve(self.K,self.R)

                for each in self.solution:
                        each[0] = round(each[0],precision)
                        if abs(each[0]) < abs(effective_zero):
                                each[0]=0.
                
                for i in range(len(self.solution)):
                        if int(i/3.) == i/3.:
                                print self.names[i][0], "... ", self.solution[i:i+3].T[0].tolist()

                
