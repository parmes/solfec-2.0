# test all Solfec-2.0 input commands

output = ARGV()
print(output)

matnum = MATERIAL (1E3, 1E9, 0.25, 0.0)
print(matnum)

nodes = [0, 0, 0,
         1, 0, 0,
	 1, 1, 0,
	 0, 1, 0,
	 0, 0, 1,
	 1, 0, 1,
	 1, 1, 1,
	 0, 1, 1]
elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]
bodnum = MESH (nodes, elements, matnum, colors)
print(bodnum)
