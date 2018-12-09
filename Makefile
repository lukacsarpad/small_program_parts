CFLAGS = -O2 -D_FILE_OFFSET_BITS=64 -DLARGEFILE_SOURCE -D_GNU_SOURCE
#link for c++, do not include electric fence
LDLIBS = -lm -lstdc++ -lboost_program_options
# link for c++
INC = 


#CC = clang $(CFLAGS) $(INC)
#CXX = clang++ -std=c++11 $(CFLAGS) $(INC)
CC = gcc $(CFLAGS) $(INC)
CXX = g++ -std=c++11 $(CFLAGS) $(INC)

LDFLAGS = -s
#LD = g++ $(LDFLAGS)
#LDFLAGS = $(LIBS) -ggdb3

listprogs=procinp_example
all:	$(listprogs)

clean:
	rm -f $(listprogs)
	rm -f *.o *~

distclean:
	rm -f *.o *~

procinp_example:	procinp_example.o procinp.o

