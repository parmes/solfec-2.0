# Solfec-2.0 input command test: SPLINE
import sys, os

d0 = os.path.dirname(os.path.realpath(sys.argv[1]))
spl0 = SPLINE (os.path.join(d0,'spline.txt'));
spl1 = SPLINE (os.path.join(d0,'spline.txt'), cache = 10)

lst2 = [0, 10, 1, 11, 2, 12, 3, 13, 4, 14, 5, 15, 6, 16];
spl2 = SPLINE (lst2);

lst3 = [[0, 10], [1, 11], [2, 12], [3, 13], [4, 14], [5, 15], [6, 16]];
spl3 = SPLINE (lst3);

print_SPLINE(spl0)
print_SPLINE(spl1)
print_SPLINE(spl2)
print_SPLINE(spl3)
