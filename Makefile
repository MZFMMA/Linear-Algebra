all: final

final: main.o data_input.o jora.o
	g++ main.o data_input.o jora.o -o final -O3

main.o: main.cpp
	g++ -c main.cpp -O3

data_input.o: data_input.cpp
	g++ -c data_input.cpp -O3

jora.o: jora.cpp
	g++ -c jora.cpp -O3

clean:
	rm -rf *.o final