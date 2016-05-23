all: bin/main.o bin/info_work.o bin/modified.o bin/segments.o bin/cysproc.o bin/keymap.o
	g++ -std=c++11 -o prog bin/main.o bin/info_work.o bin/modified.o bin/segments.o bin/cysproc.o bin/keymap.o
bin/main.o: cpp/main.cpp headers/headers.h
	g++ -std=c++11 -c cpp/main.cpp -o bin/main.o
bin/info_work.o: cpp/info_work.cpp headers/headers.h
	g++ -std=c++11 -c cpp/info_work.cpp -o bin/info_work.o
bin/modified.o: cpp/modified.cpp headers/headers.h
	g++ -std=c++11 -c cpp/modified.cpp -o bin/modified.o
bin/segments.o: cpp/segments.cpp headers/headers.h
	g++ -std=c++11 -c cpp/segments.cpp -o bin/segments.o
bin/cysproc.o: cpp/cysproc.cpp headers/headers.h
	g++ -std=c++11 -c cpp/cysproc.cpp -o bin/cysproc.o
bin/keymap.o: cpp/keymap.cpp headers/headers.h
	g++ -std=c++11 -c cpp/keymap.cpp -o bin/keymap.o
clean:
	rm *.o -f

$(shell mkdir -p bin)