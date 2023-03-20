all: final

final: main.o data_input.o ev.o
	g++ main.o data_input.o ev.o -o final -lpthread -g

main.o: main.cpp
	g++ -c main.cpp -pthread -g

data_input.o: data_input.cpp
	g++ -c data_input.cpp -lpthread -g

jora.o: ev.cpp
	g++ -c ev.cpp -pthread -g 

clean:
	rm -rf *.o final