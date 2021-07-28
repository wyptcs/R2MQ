CC = g++
SRCS = auxiliary.cpp topk_third.cpp spacetest.cpp
OBJS = $(SRCS:.cpp=.o)
TEST = test
CPPFLAGS=-std=c++11 -O3 -DNDEBUG -march=native
LIBS=-lsdsl -ldivsufsort -ldivsufsort64
LIB_DIRS = -L ~/lib
INC = -I ~/include

$(TEST) :
	$(CC) $(CPPFLAGS) $(SRCS) -o spacetest $(INC) $(LIB_DIRS) $(LIBS)

clean :
	rm -f $(TEST)
