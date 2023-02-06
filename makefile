OS := $(shell uname)
ifeq ($(OS),Darwin)
		CC      = /usr/local/bin/g++-9
        CFLAGS  = -O3 -mavx -std=c++14 -w -march=native -D $(distFunc)
        LDFLAGS = 
else
        CC      = g++
        CFLAGS  = -O3 -mavx -std=c++14 -w -march=native -D $(distFunc)
        LDFLAGS = 
endif


ifeq ($(distFunc), L2)
	SOURCES = tree/AVtree.cpp
else
	SOURCES = tree/AVtreeString.cpp
endif
OBJECTS = $(SOURCES:.cpp=.o)
	

all: main_range main_knn

main_knn: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) main_knn.cpp -o knn $(LDADD)

main_range: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) main_range.cpp -o range $(LDADD)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm -rf tree/*.o
	@rm -rf knn
	@rm -rf range
