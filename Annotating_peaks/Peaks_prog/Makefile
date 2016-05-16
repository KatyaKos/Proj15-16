all: main.o info_work.o modified.o segments.o cysproc.o
	g++ -std=c++11 -o prog main.o info_work.o modified.o segments.o cysproc.o
main.o: main.cpp headers.h
	g++ -std=c++11 -c main.cpp -o main.o
info_work.o: info_work.cpp headers.h
	g++ -std=c++11 -c info_work.cpp -o info_work.o
modified.o: modified.cpp headers.h
	g++ -std=c++11 -c modified.cpp -o modified.o
segments.o: segments.cpp headers.h
	g++ -std=c++11 -c segments.cpp -o segments.o
cysproc.o: cysproc.cpp headers.h
	g++ -std=c++11 -c cysproc.cpp -o cysproc.o
clean:
	rm *.o -f