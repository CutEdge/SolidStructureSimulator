from Solids import *

"""
Example:
          |F3     |F5     |F7
  <--a--->|       |       |
          \/      \/      \/
n1|---|---|---|---|---|---|---|---|n9
   \    l2  l3      l5  l6  l7     \ l9
  l0\                               \ (hinge axis in Z)
     (n0)                           (n10)

This triple-loaded beam with hinges matches the ASME excel calculator for a simply-supported beam at any rotation.
"""

E = 29.0e6
G = 11.3e6
I = 1./12*6*.5**3
J = 2*I
A = 6*.5
a = 2
f = -100
ax = array([-1,5,0])
ang = 0
z1=.1


nodes = []
nodes.append(Node(i=0,r=rotate([0,0,z1],ax,ang),d=[0,0,0],t=[0,0,0]))
nodes.append(Node(i=1,r=rotate([0,0,0],ax,ang)))
nodes.append(Node(i=2,r=rotate([1,0,0],ax,ang)))
nodes.append(Node(i=3,r=rotate([2,0,0],ax,ang),F=rotate([0,f,0],ax,ang)))
nodes.append(Node(i=4,r=rotate([3,0,0],ax,ang)))
nodes.append(Node(i=5,r=rotate([4,0,0],ax,ang),F=rotate([0,f,0],ax,ang)))
nodes.append(Node(i=6,r=rotate([5,0,0],ax,ang)))
nodes.append(Node(i=7,r=rotate([6,0,0],ax,ang),F=rotate([0,f,0],ax,ang)))
nodes.append(Node(i=8,r=rotate([7,0,0],ax,ang)))
nodes.append(Node(i=9,r=rotate([8,0,0],ax,ang)))
nodes.append(Node(i=10,r=rotate([8,0,z1],ax,ang),d=[0,0,0],t=[0,0,0]))

links = []
for l in range(10):
        links.append(Link(l,l,l+1,A,I,J,E,G))
links[0].I = 100
links[0].A = 100
links[0].J = 200
links[9].I = 100
links[9].A = 100
links[9].J = 200

S = System(nodes,links,verbose=True)

for l in range(10):
        if l==0:
                S.connections.append(Connection(node=0, link=0, xtype="tx", project="all"))
                S.connections.append(Connection(node=0, link=0, xtype="dx", project="all"))
                S.connections.append(Connection(node=1, link=0, xtype="tx", project="transverse"))
                S.connections.append(Connection(node=1, link=0, xtype="dx", project="all"))
        elif l==9:
                S.connections.append(Connection(node=9, link=9, xtype="tx", project="transverse"))
                S.connections.append(Connection(node=9, link=9, xtype="dx", project="all"))
                S.connections.append(Connection(node=10, link=9, xtype="tx", project="all"))
                S.connections.append(Connection(node=10, link=9, xtype="dx", project="all"))
        else:
                S.connections.append(Connection(node=l, link=l, xtype="tx", project="all"))
                S.connections.append(Connection(node=l, link=l, xtype="dx", project="all"))
                S.connections.append(Connection(node=l+1, link=l, xtype="tx", project="all"))
                S.connections.append(Connection(node=l+1, link=l, xtype="dx", project="all"))



S.buildMatrix()
saveMatrix(S.K,"matrix.txt","\t")
solution = solve(S.K,S.R)
for i in range(len(solution)):
        if int(i/3.) == i/3.:
                print S.names[i][0], "... ", rotate(solution[i:i+3].T,-ax,ang)

actual = -9.93103e-4
print  "Error: ", ( actual - rotate(solution[-3*(10-a):-3*(10-a)+3].T,-ax,ang)[0][1] ) / actual


