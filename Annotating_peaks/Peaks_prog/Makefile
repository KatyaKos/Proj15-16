all: bin/main.o bin/printer.o bin/modchain.o bin/nonmodChain.o bin/antibody.o bin/chain.o bin/reader.o
	g++ -std=c++11 -o prog bin/main.o bin/printer.o bin/modchain.o bin/nonmodChain.o bin/reader.o bin/antibody.o bin/chain.o
bin/main.o: cpp/main.cpp headers/headers.h
	g++ -std=c++11 -c cpp/main.cpp -o bin/main.o
bin/printer.o: cpp/printer.cpp headers/headers.h
	g++ -std=c++11 -c cpp/printer.cpp -o bin/printer.o
bin/reader.o: cpp/reader.cpp headers/headers.h
	g++ -std=c++11 -c cpp/reader.cpp -o bin/reader.o
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