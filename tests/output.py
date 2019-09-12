# Solfec-2.0 input command test: OUTPUT

# define print default outputs
OUTPUT()
print_OUTPUTS()

# add bodies
matnum = MATERIAL (1E3, 1E9, 0.25, 0.0)
center = (0, 0, 0)
radius = (1, 2, 3)
bodnum0 = ELLIP (center, radius, matnum, 1)
center = (1, 1, 1)
radius = (3, 2, 1)
bodnum1 = ELLIP (center, radius, matnum, 2)

# define subset OUTPUTS
OUTPUT(['COLOR', 'DISPL'], bodnum0)
OUTPUT(['LINVEL', 'STRESS'], bodnum1, modes=['MESH'])
OUTPUT(['CF', 'AREA', 'CPAIR'], [bodnum0, bodnum1], modes='CD')

# print all outputs
print_OUTPUTS()
