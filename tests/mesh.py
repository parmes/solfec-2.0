# Solfec-2.0 input command test: MESH

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

print_MESH(bodnum)
