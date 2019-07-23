# Solfec-2.0 input command test: GRAVITY

spl = SPLINE ([0, 0, -5, 1, 2, -9.81, 3, -9.81])

def call(t): return 0.0

GRAVITY (0.0, call, spl)

print_GRAVITY()
