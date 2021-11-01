import numpy as np
import scipy
#from scipy.sparse.linalg import spsolve 
from scipy.sparse import csr_matrix, coo_matrix
import cupy as cp
from cupyx.scipy.sparse import spmatrix as csr_gpu
from cupyx.scipy.sparse.linalg import spsolve 
import time

t1 = time.time()

dim = 3

poisson = 0.3
youngModulus = 2000


class Triplet:
  def __init__(self, row,col,val):
    self.row = row
    self.col = col
    self.val = val

class Element:
  def __init__(self,node1,node2,node3,node4):
    self.nodesIds = np.zeros(4, dtype = int)
    self.nodesIds[0] = node1
    self.nodesIds[1] = node2
    self.nodesIds[2] = node3
    self.nodesIds[3] = node4
    self.B = np.zeros((6,12))

  def CalculateStiffnessMatrix(self,D):
    x = np.array([nodes_x[self.nodesIds[0]], nodes_x[self.nodesIds[1]], nodes_x[self.nodesIds[2]], nodes_x[self.nodesIds[3]]])
    y = np.array([nodes_y[self.nodesIds[0]], nodes_y[self.nodesIds[1]], nodes_y[self.nodesIds[2]], nodes_y[self.nodesIds[3]]])
    z = np.array([nodes_z[self.nodesIds[0]], nodes_z[self.nodesIds[1]], nodes_z[self.nodesIds[2]], nodes_z[self.nodesIds[3]]])
    C = np.matrix([np.ones(4),x,y,z])
    C = C.T
    IC = C.I
    for i in range(4):
      self.B[0, 3 * i + 0] = IC[1,i]
      self.B[0, 3 * i + 1] = 0.
      self.B[0, 3 * i + 2] = 0.

      self.B[1, 3 * i + 0] = 0.
      self.B[1, 3 * i + 1] = IC[2,i]
      self.B[1, 3 * i + 2] = 0.

      self.B[2, 3 * i + 0] = 0.
      self.B[2, 3 * i + 1] = 0.
      self.B[2, 3 * i + 2] = IC[3,i]

      self.B[3, 3 * i + 0] = IC[2,i]
      self.B[3, 3 * i + 1] = IC[1,i]
      self.B[3, 3 * i + 2] = 0.

      self.B[4, 3 * i + 0] = 0.
      self.B[4, 3 * i + 1] = IC[3,i]
      self.B[4, 3 * i + 2] = IC[2,i]

      self.B[5, 3 * i + 0] = IC[3,i]
      self.B[5, 3 * i + 1] = 0.
      self.B[5, 3 * i + 2] = IC[1,i]
    
    K = np.matmul(np.matmul(self.B.T, D), self.B) * (np.linalg.det(C) / 6)

    triplets = []
    for i in range(4):
      for j in range(4):
        trplt11 = Triplet(3 * self.nodesIds[i] + 0, 3 * self.nodesIds[j] + 0, K[3*i + 0, 3*j + 0])
        trplt12 = Triplet(3 * self.nodesIds[i] + 0, 3 * self.nodesIds[j] + 1, K[3*i + 0, 3*j + 1])
        trplt13 = Triplet(3 * self.nodesIds[i] + 0, 3 * self.nodesIds[j] + 2, K[3*i + 0, 3*j + 2])
        trplt21 = Triplet(3 * self.nodesIds[i] + 1, 3 * self.nodesIds[j] + 0, K[3*i + 1, 3*j + 0])
        trplt22 = Triplet(3 * self.nodesIds[i] + 1, 3 * self.nodesIds[j] + 1, K[3*i + 1, 3*j + 1])
        trplt23 = Triplet(3 * self.nodesIds[i] + 1, 3 * self.nodesIds[j] + 2, K[3*i + 1, 3*j + 2])
        trplt31 = Triplet(3 * self.nodesIds[i] + 2, 3 * self.nodesIds[j] + 0, K[3*i + 2, 3*j + 0])
        trplt32 = Triplet(3 * self.nodesIds[i] + 2, 3 * self.nodesIds[j] + 1, K[3*i + 2, 3*j + 1])
        trplt33 = Triplet(3 * self.nodesIds[i] + 2, 3 * self.nodesIds[j] + 2, K[3*i + 2, 3*j + 2])
  
        if trplt11.val != 0:       
          triplets.append(trplt11)
        if trplt12.val != 0:       
          triplets.append(trplt12)
        if trplt13.val != 0:       
          triplets.append(trplt13)
        if trplt21.val != 0:       
          triplets.append(trplt21)
        if trplt22.val != 0:       
          triplets.append(trplt22)
        if trplt23.val != 0:       
          triplets.append(trplt23)
        if trplt31.val != 0:       
          triplets.append(trplt31)
        if trplt32.val != 0:       
          triplets.append(trplt32)
        if trplt33.val != 0:       
          triplets.append(trplt33)


    return triplets

class Constraint:

    def __init__(self, node, t):

      self.node = node
      self.type = t

#Initialization   

with open('nodes.txt','r') as nodes:
  nodes_count = int(nodes.readline())
  nodes_x = np.zeros(nodes_count)
  nodes_y = np.zeros(nodes_count)
  nodes_z = np.zeros(nodes_count)
  i = 0
  for line in nodes:
    x,y,z = map(float,line.split())
    nodes_x[i] = x
    nodes_y[i] = y
    nodes_z[i] = z
    i += 1

with open('elements.txt','r') as el:

  elements_size = int(el.readline())
  elements = []
  
  for line in el:
    n1,n2,n3,n4 = map(int,line.split())
    e = Element(n1,n2,n3,n4)
    elements.append(e)
    
with open('constraints.txt','r') as cons:
  constraints_size = int(cons.readline())
  constraints = []
  for line in cons:
    node,type = map(int,line.split())
    c = Constraint(node,type)
    constraints.append(c)

with open('loads.txt','r') as l:
  loads_count = int(l.readline())
  loads = np.zeros(dim*nodes_count)
  for line in l:
    node = int(line.split()[0])
    x,y,z = map(float,line.split()[1:])
    loads[dim*node + 0] = x
    loads[dim*node + 1] = y
    loads[dim*node + 2] = z

print(f'Load data: {time.time() - t1}')




#D matrix

D = np.zeros((6,6))
D[0,0] = 1.
D[1,1] = 1.
D[2,2] = 1.
D[3,3] = (1 - 2 * poisson) / (2 * (1 - poisson))
D[4,4] = (1 - 2 * poisson) / (2 * (1 - poisson))
D[5,5] = (1 - 2 * poisson) / (2 * (1 - poisson))
D[0,1] = poisson / (1 - poisson)
D[0,2] = poisson / (1 - poisson)
D[1,0] = poisson / (1 - poisson)
D[1,2] = poisson / (1 - poisson)
D[2,0] = poisson / (1 - poisson)
D[2,1] = poisson / (1 - poisson)
D *= youngModulus * (1 - poisson) / ((1 + poisson)*(1 - 2 * poisson))

triplets = []

for element in elements:
  triplets = triplets + element.CalculateStiffnessMatrix(D)

N = nodes_count * dim

global_K_rows = []
global_K_cols = []
global_K_values = []

for triplet in triplets:
  for i,j,k in zip(global_K_rows, global_K_cols, range(len(global_K_values))):
    if triplet.row == i and triplet.col == j:
      global_K_values[k] += triplet.val
      break
  else:
    if triplet.row == 0 and triplet.col == 0:
      print(triplet.val)
    global_K_rows.append(triplet.row)
    global_K_cols.append(triplet.col)
    global_K_values.append(triplet.val)

print(f'Create global_K {time.time() - t1}')


global_K = coo_matrix((global_K_values, (global_K_rows, global_K_cols)), shape = (N,N))
#print(global_K)

def ApplyConstraints(K, constraints):
  indToConstraint = []
  for i in constraints:
    if i.type == 1:
      indToConstraint.append(3 * i.node + 0)
    if i.type == 2:
      indToConstraint.append(3 * i.node + 1)
    if i.type == 3:
      indToConstraint.append(3 * i.node + 2)
    if i.type == 4:
      indToConstraint.append(3 * i.node + 0)
      indToConstraint.append(3 * i.node + 1)
    if i.type == 5:
      indToConstraint.append(3 * i.node + 0)
      indToConstraint.append(3 * i.node + 2)
    if i.type == 6:
      indToConstraint.append(3 * i.node + 1)
      indToConstraint.append(3 * i.node + 2)
    if i.type == 7:
      indToConstraint.append(3 * i.node + 0)
      indToConstraint.append(3 * i.node + 1)
      indToConstraint.append(3 * i.node + 2)
      
    
  for i,j,k in zip(K.row,K.col,range(len(K.data))):
    for it in indToConstraint:
      if i == j and i == it:
        K.data[k] = 1
        continue
      if i == it:
        K.data[k] = 0
      if j == it:
        K.data[k] = 0
  K_new = csr_matrix((K.data, (K.row, K.col)), shape = K.shape)
  return K_new


global_K = ApplyConstraints(global_K,constraints)

global_K_gpu = csr_gpu(global_K)
loads_gpu = cp.array(loads)
x_gpu = spsolve(global_K_gpu, loads_gpu)
print(x_gpu)

print(f'Apply constraints: {time.time() - t1}')

print(spsolve(global_K,loads))

print(f'Solve {time.time() - t1}')

t2 = time.time()

print(f'Time: {t2 - t1}')
    

