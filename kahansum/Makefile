#CFLAGS = -msse2 -mfpmath=sse,387 -O2 -march=nocona\
#	    -mtune=nocona # -ggdb3
CFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -DLARGEFILE_SOURCE -D_GNU_SOURCE -g
#link for c++, do not include electric fence
LDLIBS = -lm -lstdc++
# link for c++, also use a malloc debugger for testing purposes
#LDLIBS = libsdf.a -lm -lstdc++ -lefence
INC = 


CC = g++ $(CFLAGS) $(INC)
CXX = g++ -std=c++17 $(CFLAGS) $(INC)

LDFLAGS =
#LD = g++ $(LDFLAGS)
#LDFLAGS = $(LIBS) -ggdb3

probalist = $(basename $(wildcard test*.cc))
all:	$(probalist)

clean:
	rm -f $(probalist)
	rm -f *.o *~ *.h5

distclean:
	rm -f *.o *~



testkahansum:	testkahansum.cc

