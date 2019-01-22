all: InterLig

InterLig: InterLig.cpp InterLig.h

	g++ -static -O3 -ffast-math -lm -o InterLig InterLig.cpp 
clean:
	rm -f InterLig

test:
	sh ./test.sh
