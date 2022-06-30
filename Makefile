CC=mpic++
CFLAGS=-c -O3
LDFLAGS=
SOURCES=main2.cpp solve2.cpp norma.cpp ent.cpp outmatrix.cpp
OBJECTS=$(SOURCES:.cpp=.o)
DEPS=$(SOURSES:.cpp=.d)
EXECUTABLE=output

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS)

