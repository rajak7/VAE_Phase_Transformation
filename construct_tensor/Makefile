c_feature: createfeature.o atom_property.o
	gcc -o c_feature createfeature.o atom_property.o -lm

createfeature.0: createfeature.c atom_property.h createfeature.h
	gcc  -c -Wall -std=c99 createfeature.c

atom_property.o: atom_property.c atom_property.h
	gcc  -c -Wall -std=c99 atom_property.c

clean:
	rm -f *.o c_feature
