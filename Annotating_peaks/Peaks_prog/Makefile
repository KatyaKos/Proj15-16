all: main.o info_work.o modified.o segments.o cysproc.o keymap.o
	g++ -std=c++11 -o prog main.o info_work.o modified.o segments.o cysproc.o keymap.o
main.o: main.cpp headers/headers.h
	g++ -std=c++11 -c main.cpp -o main.o
info_work.o: info_work.cpp headers/headers.h
	g++ -std=c++11 -c info_work.cpp -o info_work.o
modified.o: modified.cpp headers/headers.h
	g++ -std=c++11 -c modified.cpp -o modified.o
segments.o: segments.cpp headers/headers.h
	g++ -std=c++11 -c segments.cpp -o segments.o
cysproc.o: cysproc.cpp headers/headers.h
	g++ -std=c++11 -c cysproc.cpp -o cysproc.o
keymap.o: keymap.cpp headers/headers.h
	g++ -std=c++11 -c keymap.cpp -o keymap.o
clean:
	rm *.o -f