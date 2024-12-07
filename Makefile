execute: hamming.o
	g++ hamming.o -o execute

hamming.o: Hamming.cpp Hamming.h
	g++ -c Hamming.cpp -o hamming.o

run: execute
	./execute

clean:
	rm -f *.o execute
