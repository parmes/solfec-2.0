# Solfec-2.0 python routine --> MESH_HEX: create a tri-linear hexahedral mesh based particle

def SUB (a, b, c):
  c[0] = a[0] - b[0]
  c[1] = a[1] - b[1]
  c[2] = a[2] - b[2]

def PRODUCT (a, b, c):
  c[0] = a[1]*b[2] - a[2]*b[1]
  c[1] = a[2]*b[0] - a[0]*b[2]
  c[2] = a[0]*b[1] - a[1]*b[0]

def DOT (a, b):
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def COPY (a, b):
  b[0] = a[0]
  b[1] = a[1]
  b[2] = a[2]

# linear shape functions for the hexahedron
def HEX0(x,y,z): return (0.125*(1.0-(x))*(1.0-(y))*(1.0-(z)))
def HEX1(x,y,z): return (0.125*(1.0+(x))*(1.0-(y))*(1.0-(z)))
def HEX2(x,y,z): return (0.125*(1.0+(x))*(1.0+(y))*(1.0-(z)))
def HEX3(x,y,z): return (0.125*(1.0-(x))*(1.0+(y))*(1.0-(z)))
def HEX4(x,y,z): return (0.125*(1.0-(x))*(1.0-(y))*(1.0+(z)))
def HEX5(x,y,z): return (0.125*(1.0+(x))*(1.0-(y))*(1.0+(z)))
def HEX6(x,y,z): return (0.125*(1.0+(x))*(1.0+(y))*(1.0+(z)))
def HEX7(x,y,z): return (0.125*(1.0-(x))*(1.0+(y))*(1.0+(z)))

# nodes := [(x0,y0,z10), (x1,y1,z1), ..., (x7,y7,z7)] defininf eight corners
# i, j, k := integer numbers of subdivisions along x, y, and z
# material := material number
# colors := [c0, c1, ..., c7] list of face colors
# dx := [x0, x1, ..., xi] where x0,x1,... > 0 and SUM{x0+x1+...} = 1.0, define ratios of i subdivisions along x
# dy, dz := are defined in analogy to dx
def MESH_HEX (nodes, i, j, k, material, colors, dx = None, dy = None, dz = None):
  from solfec import MESH
  vxvy = [0,0,0];
  vx = [0,0,0];
  vy = [0,0,0];
  vz = [0,0,0];

  SUB (nodes [1], nodes [0], vx);
  SUB (nodes [2], nodes [1], vy);
  SUB (nodes [4], nodes [0], vz);
  PRODUCT (vx, vy, vxvy);
  if (DOT (vxvy, vz) < 0):  # change orientation: swap upper and lowe nodes
    copy = [[0,0,0]]*4;
    for ii in range(0,4): COPY (nodes [ii], copy [ii]);
    for ii in range(0,4): COPY (nodes [ii+4], nodes [ii]);
    for ii in range(0,4): COPY (copy [ii], nodes [ii+4]);

  nod = [0]*(3*(i+1)*(j+1)*(k+1));
  ele = [0]*(10*i*j*k);
  col = [0]*(1+6*((2*i*j)+(2*j*k)+(2*i*k)));

  # create the unit cube mesh
  ddx =  2.0 / float(i);
  ddy =  2.0 / float(j);
  ddz =  2.0 / float(k);
  mx = my = mz = 0;

  if (dx == None):
    dx = [ddx]*i;
  else:
    x = 0.0;
    for ii in range (0, i): x += dx[ii];
    for ii in range (0, i): dx[ii] *= 2.0/x;

  if (dy == None):
    dy = [ddy]*j;
  else:
    y = 0.0;
    for jj in range (0, j): y += dy[jj];
    for jj in range (0, j): dy[jj] *= 2.0/y;

  if (dz == None):
    dz = [ddz]*k;
  else:
    z = 0.0;
    for kk in range (0, k): z += dz[kk];
    for kk in range (0, k): dz[kk] *= 2.0/z;

  # nodes
  z = -1.0;
  nn = 0;
  for kk in range(0, k+1): 
    y = -1.0;
    for jj in range(0, j+1):
      x = -1.0;
      for ii in range(0, i+1):
        nod [3*nn+0] = nodes[0][0]*HEX0(x,y,z) + nodes[1][0]*HEX1(x,y,z) + nodes[2][0]*HEX2(x,y,z) + nodes[3][0]*HEX3(x,y,z)+\
                       nodes[4][0]*HEX4(x,y,z) + nodes[5][0]*HEX5(x,y,z) + nodes[6][0]*HEX6(x,y,z) + nodes[7][0]*HEX7(x,y,z);
        nod [3*nn+1] = nodes[0][1]*HEX0(x,y,z) + nodes[1][1]*HEX1(x,y,z) + nodes[2][1]*HEX2(x,y,z) + nodes[3][1]*HEX3(x,y,z)+\
                       nodes[4][1]*HEX4(x,y,z) + nodes[5][1]*HEX5(x,y,z) + nodes[6][1]*HEX6(x,y,z) + nodes[7][1]*HEX7(x,y,z);
        nod [3*nn+2] = nodes[0][2]*HEX0(x,y,z) + nodes[1][2]*HEX1(x,y,z) + nodes[2][2]*HEX2(x,y,z) + nodes[3][2]*HEX3(x,y,z)+\
                       nodes[4][2]*HEX4(x,y,z) + nodes[5][2]*HEX5(x,y,z) + nodes[6][2]*HEX6(x,y,z) + nodes[7][2]*HEX7(x,y,z);
        x += dx[min(ii,i-1)];
        nn += 1;
      y += dy[min(jj,j-1)];
    z += dz[min(kk,k-1)];

  # elements
  face = [[1,4,3,2], [1,2,6,5], [2,3,7,6], [3,4,8,7], [1,5,8,4], [5,6,7,8]];
  n = (i+1)*(j+1); 
  col[0] = colors[0]; # will not be used but let it be set sensibly
  ee = 0;
  cc = 1;
  for kk in range (0, k):
    for jj in range (0, j):
      for ii in range (0, i):
        nn = ((kk*n)+(jj*(i+1))+ii);

        ele[ee+0] = 8;
        ele[ee+1] = nn;
        ele[ee+2] = nn+1;
        ele[ee+3] = nn+i+2;
        ele[ee+4] = nn+i+1;
        ele[ee+5] = ele[ee+1]+n;
        ele[ee+6] = ele[ee+2]+n;
        ele[ee+7] = ele[ee+3]+n;
        ele[ee+8] = ele[ee+4]+n;
        ele[ee+9] = material;

        if (kk == 0): # face 1
          col[cc+0] = 4;
          for m in range (0, 4): col[cc+m+1] = ele[ee+face[0][m]];
          col[cc+5] = colors[0];
          cc += 6;
        if (jj == 0): # face 2
          col[cc+0] = 4;
          for m in range(0, 4): col[cc+m+1] = ele[ee+face[1][m]];
          col[cc+5] = colors[4];
          cc += 6;
        if (ii == (i-1)): # face 3
          col[cc+0] = 4;
          for m in range(0, 4): col[cc+m+1] = ele[ee+face[2][m]];
          col[cc+5] = colors[3];
          cc += 6;
        if (jj == (j-1)): # face 4
          col[cc+0] = 4;
          for m in range(0, 4): col[cc+m+1] = ele[ee+face[3][m]];
          col[cc+5] = colors[2];
          cc += 6;
        if (ii == 0): # face 5
          col[cc+0] = 4;
          for m in range(0, 4): col[cc+m+1] = ele[ee+face[4][m]];
          col[cc+5] = colors[1];
          cc += 6;
        if (kk == (k-1)): # face 6
          col[cc+0] = 4;
          for m in range(0, 4): col[cc+m+1] = ele[ee+face[5][m]];
          col[cc+5] = colors[5];
          cc += 6;

        ee += 10

  return MESH (nod, ele, material, col)
