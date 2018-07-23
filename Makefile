CXX = mpiCC
CXXFLAGS = -std=c++0x -Wall -Wextra -Werror -Wshadow -O3 -DNDEBUG

INCLUDES =
LDFLAGS =
LIBS =

# blas
INCLUDES += -I/usr/lib64/atlas/include/
LDFLAGS += -L/usr/lib64/atlas/
LIBS += -lcblas -latlas

# likwid
CXXFLAGS += -DUSE_LIKWID -pthread
INCLUDES += -I/usr/local/likwid/include/
LDFLAGS += -L/usr/local/likwid/lib/
LIBS += -llikwid

TARGET = matmult
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

$(TARGET).o: $(TARGET).cpp matrixmult.cpp matrixmult.h Timer.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(TARGET).cpp

clean:
	@$(RM) -rf *.o $(TARGET)
