#CFLAGS = -msse2 -mfpmath=sse,387 -O2 -march=nocona\
#	    -mtune=nocona # -ggdb3
CFLAGS = -Wall -O2 -D_FILE_OFFSET_BITS=64 -DLARGEFILE_SOURCE
# -g
#link for c++, do not include electric fence
LDLIBS = -lm -lstdc++
# link for c++, also use a malloc debugger for testing purposes
#LDLIBS = libsdf.a -lm -lstdc++ -lefence
INC = 


CC = h5c++ $(CFLAGS) $(INC)
CXX = h5c++ -std=c++11 $(CFLAGS) $(INC)

LDFLAGS =
#LD = g++ $(LDFLAGS)
#LDFLAGS = $(LIBS) -ggdb3

probalist = $(basename $(wildcard proba*.cc))
all:	$(probalist)

clean:
	rm -f $(probalist)
	rm -f *.o *~ *.h5

distclean:
	rm -f *.o *~



proba1:		proba1.o alhdf5.o 

proba2:		proba2.o alhdf5.o 

proba3:		proba3.o alhdf5.o 

