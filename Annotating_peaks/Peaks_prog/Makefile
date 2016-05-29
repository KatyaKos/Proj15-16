all: bin/main.o bin/modchain.o bin/nonmodChain.o bin/antibody.o bin/chain.o bin/spectrum.o
	g++ -std=c++11 -o prog bin/main.o bin/modchain.o bin/nonmodChain.o bin/spectrum.o bin/antibody.o bin/chain.o
bin/main.o: cpp/main.cpp headers/headers.h
	g++ -std=c++11 -c cpp/main.cpp -o bin/main.o
bin/spectrum.o: cpp/spectrum.cpp headers/headers.h
	g++ -std=c++11 -c cpp/spectrum.cpp -o bin/spectrum.o
bin/antibody.o: cpp/antibody.cpp headers/headers.h
	g++ -std=c++11 -c cpp/antibody.cpp -o bin/antibody.o
bin/chain.o: cpp/chain.cpp headers/headers.h
	g++ -std=c++11 -c cpp/chain.cpp -o bin/chain.o
bin/modchain.o: cpp/modchain.cpp headers/headers.h
	g++ -std=c++11 -c cpp/modchain.cpp -o bin/modchain.o
bin/nonmodChain.o: cpp/nonmodChain.cpp headers/headers.h
	g++ -std=c++11 -c cpp/nonmodChain.cpp -o bin/nonmodChain.o
clean:
	rm bin/*.o -f

$(shell mkdir -p bin)