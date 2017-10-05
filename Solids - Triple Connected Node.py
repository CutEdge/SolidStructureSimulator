from Solids import *

"""
Example:
          |F
          |
    l0    \/
n0|---|---|n2
     /  l1
    /l2
   /
   n3
"""

E = 29.0e6
G = 11.3e6
I = 1./12*6*.5**3
J = 2*I
A = 6*.5

nodes = []
nodes.append(Node(0,r=[0,0,0],t=[0,0,None],d=[0,0,0]))
nodes.append(Node(1,r=[5,0,0]))
nodes.append(Node(2,r=[10,0,0],F=[0,-100,0]))
nodes.append(Node(3,r=[0,-5,0],t=[0,0,None],d=[0,0,0]))

links = []
links.append(Link(0,0,1,A,I,J,E,G))
links.append(Link(1,1,2,A,I,J,E,G))
links.append(Link(2,3,1,A,I,J,E,G))

S = System(nodes,links)

S.connections.append(Connection(0,0, xtype="tx", project="all"))
S.connections.append(Connection(0,0, xtype="dx", project="all"))

S.connections.append(Connection(1,0, xtype="tx", project="all"))
S.connections.append(Connection(1,0, xtype="dx", project="all"))
S.connections.append(Connection(1,1, xtype="tx", project="all"))
S.connections.append(Connection(1,1, xtype="dx", project="all"))
S.connections.append(Connection(1,2, xtype="tx", project="all"))
S.connections.append(Connection(1,2, xtype="dx", project="all"))

S.connections.append(Connection(2,1, xtype="tx", project="all"))
S.connections.append(Connection(2,1, xtype="dx", project="all"))

S.connections.append(Connection(3,2, xtype="tx", project="all"))
S.connections.append(Connection(3,2, xtype="dx", project="all"))


S.buildMatrix()
#saveMatrix(S.K,"matrix.txt","\t")
solution = solve(S.K,S.R)
for i in range(len(solution)):
        if int(i/3.) == i/3.:
                print S.names[i][0], "... ", solution[i:i+3].T
