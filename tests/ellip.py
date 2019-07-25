# Solfec-2.0 input command test: ELLIP

matnum = MATERIAL (1E3, 1E9, 0.25, 0.0)

center = (0, 0, 0)
radius = (1, 2, 3)
bodnum0 = ELLIP (center, radius, matnum, 1)

center = (1, 1, 1)
radius = (3, 2, 1)
bodnum1 = ELLIP (center, radius, matnum, 2)

print_ELLIP(bodnum0)
print_ELLIP(bodnum1)
