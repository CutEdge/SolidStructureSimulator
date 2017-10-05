from Solids import *

"""
Example:
           |F
  <--a---> | <---------b---------->
          \/
n0|---|---|---|---|---|---|---|---|n8
        l1  l2      l4  l5      l7


Compare the results of this to the equation { d = F*a**2b**2 / (3*E*I*(a+b)) }
When t0 and t8 constraints are removed, a singular matrix occurs *only when beam is aligned with x axis*
Okay, the change happens when rotating it because you can't rotate the node constraints. t=[0,0,None] stays in the original coordinates.
This is why the cantilever works but the simply supported doesn't!
When you rotate a cantilever, the end conditions are the same, because all 3 coordinates are the same anyway.
But when you rotate a hinge-ended beam, the hinges don't rotate. So you're hinging about a weird axis.
But there's still a problem. Removing the t0 and t8 constraints gives error when rotated - different levels of error for different rotations.
So it's important to have the problem fully constrained, even with respect to irrelevant axes.


"""

E = 29.0e6
G = 11.3e6
I = 1./12*6*.5**3
J = 2*I
A = 6*.5
a = 2
f = -100
ax = array([-1,1,0])
ang = 0


nodes = []
for n in range(9):
        nodes.append(Node(i=n,r=rotate([n,0,0],ax,ang)))
nodes[0] = Node(i=0,r=[0,0,0],d=[0,0,0],t=[0,0,None])
nodes[a] = Node(i=a,r=rotate([a,0,0],ax,ang),F=rotate([0,f,0],ax,ang))
nodes[8] = Node(i=8,r=rotate([8,0,0],ax,ang),d=[0,0,0],t=[0,0,None])

links = []
for l in range(8):
        links.append(Link(l,l,l+1,A,I,J,E,G))

S = System(nodes,links,verbose=True)

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
                print S.names[i][0], "... ", rotate(solution[i:i+3].T,-ax,ang)

actual = f*(8.-a)**2*a**2/(3*E*I*8)
print  "Error: ", ( actual - rotate(solution[-3*(9-a):-3*(9-a)+3].T,-ax,ang)[0][1] ) / actual


