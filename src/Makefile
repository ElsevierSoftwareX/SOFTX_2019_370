CC	  = gcc
CXX	  = g++
XFLAGS    = -Wall
UFLAGS    = # User overridable at command line
STD	  = -std=c++11

INC1	  = /vlsci/VR0345/shared/pwiz 
INC2	  = /vlsci/VR0345/shared/pwiz/libraries/boost-aux
INC3	  = /usr/local/boost/1.57.0-gcc/include
INCDIRS   = -I${INC1} -I${INC2} -I${INC3}

LIBDIR1   = /usr/local/pwiz/3.0.7069-gcc/lib
LIBDIR2   = /usr/local/boost/1.57.0-gcc/lib
LIBPWIZ   = -lpwiz
LIBCHRONO = -lboost_chrono
LDFLAGS   = -L${LIBDIR1} -L${LIBDIR2} ${LIBPWIZ} ${LIBCHRONO}

CXXFLAGS  = ${STD} ${INCDIRS} ${XFLAGS} ${UFLAGS}
RM	  = rm -f

SRCS=hitime.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: hitime

hitime.out: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o hitime.out $(OBJS)

hitime.o: hitime.cpp

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) hitime.out
