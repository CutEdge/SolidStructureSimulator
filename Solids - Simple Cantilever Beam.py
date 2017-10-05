from Solids import *

"""
Example:
                                  |F
                                  |
                                  \/
n0|---|---|---|---|---|---|---|---|n8
        l1  l2      l4  l5      l7


Compare the results of this to the equation { d = F*L^3 / (3*E*I) }
"""

E = 29.0e6
G = 11.3e6
I = 1./12*6*.5**3
J = 2*I
A = 6*.5
ax = array([-1,1,0])
ang = 0


nodes = []
for n in range(9):
        nodes.append(Node(i=n,r=rotate([n,0,0],ax,ang)))
nodes[0] = Node(i=0,r=[0,0,0],d=[0,0,0],t=[0,0,0])
nodes[8] = Node(i=8,r=rotate([8,0,0],ax,ang),F=rotate([0,-100,0],ax,ang))

links = []
for l in range(8):
        links.append(Link(l,l,l+1,A,I,J,E,G))

S = System(nodes,links)

for l in range(8):
        S.connections.append(Connection(node=l, link=l, xtype="tx", project="all"))
        S.connections.append(Connection(node=l, link=l, xtype="dx", project="all"))
        S.connections.append(Connection(node=l+1, link=l, xtype="tx", project="all"))
        S.connections.append(Connection(node=l+1, link=l, xtype="dx", project="all"))



S.buildMatrix()
saveMatrix(S.K,"matrix.txt","\t")
solution = solve(S.K,S.R)
for i in range(len(solution)):
        if int(i/3.) == i/3.:
                print S.names[i][0], "... ", solution[i:i+3].T

actual = -100.*8.**3/(3*E*I)
print  "Error: ", ( actual - rotate(solution[-3:].T,-ax,ang)[0][1] ) / actual


