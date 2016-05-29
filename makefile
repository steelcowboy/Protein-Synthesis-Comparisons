CC=g++
TARGET=protein_synthesis
SRCDIR=src
CMPDIR=comparison

all: bin 
	$(CC) $(SRCDIR)/$(TARGET).cc -o bin/$(TARGET)
	cp $(SRCDIR)/$(TARGET).py bin/

comparison: bin
	$(CC) $(SRCDIR)/$(CMPDIR)/$(TARGET)_multithread.cc -lpthread -o bin/$(TARGET)_multithread
	$(CC) $(SRCDIR)/$(CMPDIR)/$(TARGET).cc -o bin/$(TARGET)
	cp $(SRCDIR)/$(CMPDIR)/$(TARGET).py bin/

bin:
	mkdir $@

clean:
	$(RM) bin/*
