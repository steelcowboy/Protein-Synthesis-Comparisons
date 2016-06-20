CC=g++

TARGET=protein_synthesis

all: bin 
	$(CC) $(TARGET)_multithread.cc -lpthread -o bin/$(TARGET)_multithread
	$(CC) $(TARGET).cc -o bin/$(TARGET)
	cp $(TARGET).py bin/

bin:
	mkdir $@

clean:
	$(RM) bin/*
