CC = g++
SRCS = DFUDSIndexing.cpp FischerHeunIndexing.cpp DFUDSEncoding.cpp BPIndexing.cpp
SRCS_DCCDFUDS = dcctest.cpp rmqFischerDFUDS-master/Basicrmq.cpp rmqFischerDFUDS-master/DFUDSrmq.cpp
SRCS_DCCBP = dccbptest.cpp rmq-master/RMQRMM64.cpp rmq-master/includes/Basic_rmq.cpp
OBJS = $(SRCS:.cpp=.o)
DCCDFUDS = dccdfuds
DCCBP = dccbp
SEA = sea
EXAMPLE = example
TUTORIAL = tutorial
CPPFLAGS=-std=c++11 -O3 -DNDEBUG -march=native
LIBS=-lsdsl -ldivsufsort -ldivsufsort64
LIB_DIRS = -L ~/lib
INC = -I ~/include
all : obj static

$(DCCDFUDS) :
	$(CC) $(CPPFLAGS) $(SRCS_DCCDFUDS) -o dcctest

$(DCCBP) :
	$(CC) $(CPPFLAGS) $(SRCS_DCCBP) -o dccbptest

$(SEA) :
	$(CC) $(CPPFLAGS) SEATEST.cpp -o SEATEST $(INC) $(LIB_DIRS) $(LIBS)

$(EXAMPLE) :
	$(CC) $(CPPFLAGS) example.cpp -o $(EXAMPLE) DFUDSR2MQ.a $(INC) $(LIB_DIRS) $(LIBS)

$(TUTORIAL) :
	$(CC) $(CPPFLAGS) tutorial.cpp -o tutorial DFUDSR2MQ.a $(INC) $(LIB_DIRS) $(LIBS)

static : $(OBJS)
	ar rc DFUDSR2MQ.a $(OBJS)

obj :
	$(CC) $(CPPFLAGS)  -c $(SRCS) $(INC) $(LIB_DIRS) $(LIBS)

clean :
	rm -f $(SEA)
	rm -f $(EXAMPLE)
	rm -f $(TUTORIAL)
	rm -f *.o
	rm -f *.a

