CC = g++
CFLAGS = -std=c++20 -Wall -g

SRCS = src/main.cpp 
HDRS = helper.hpp met.hpp

OBJS = $(SRCS:.cpp=.o) 

TARGET = program

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

.PHONY: clean

clean:
	rm -f $(TARGET) $(OBJS)

src/%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@




