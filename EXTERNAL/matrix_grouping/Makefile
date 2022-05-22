all : lib c_funcs

c_funcs: src/C/def.h src/C/dmatrix.h
	R CMD SHLIB src/C/main.cpp src/C/def.c src/C/dmatrix.c -o lib/c_funcs.so
	
clean: 
	rm -f src/C/*.o src/C/*.so lib/*.so
	
lib: 
	mkdir lib
