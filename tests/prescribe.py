# Solfec-2.0 input command test: PRESCRIBE
import sys, os

d0 = os.path.dirname(os.path.realpath(sys.argv[1]))
spl0 = SPLINE (os.path.join(d0,'spline.txt'));
lst1 = [0, 10, 1, 11, 2, 12, 3, 13, 4, 14, 5, 15, 6, 16];
spl1 = SPLINE (lst1);

matnum = MATERIAL (1E3, 1E9, 0.25, 0.0)

nodes = [0, 0, 0,
         1, 0, 0,
	 1, 1, 0,
	 0, 1, 0,
	 0, 0, 1,
	 1, 0, 1,
	 1, 1, 1,
	 0, 1, 1,
	 .5, .5, 2,
	 1, 1, 3,
	 0, 1, 3,
	 .5, .5, 3,
	 .5, .5, 4]
elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum,
            5, 4, 5, 6, 7, 8, matnum,
            6, 6, 7, 8, 9, 10, 11, matnum,
	    4, 9, 10, 11, 12, matnum]
colors = [1, 4, 0, 1, 2, 3, 2,
             4, 4, 5, 6, 7, 3,
             3, 10, 11, 12, 4]
bodnum = MESH (nodes, elements, matnum, colors)

def call1(t):
  return (1, 0, 0)

res0 = PRESCRIBE (bodnum, (1, 0, 0), linear=(0.0, 1.0, spl0))
res1 = PRESCRIBE (bodnum, [(0, 0, 0), (1, 0, 0)], linear=call1)
res2 = PRESCRIBE (bodnum, color=4, angular=(1.0, spl1, 0.0))
res3 = PRESCRIBE (bodnum, color=4, angular=call1)

print_PRESCRIBE(res0)
print_PRESCRIBE(res1)
print_PRESCRIBE(res2)
print_PRESCRIBE(res3)
