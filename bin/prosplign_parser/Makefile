CC=g++
CFLAGS=-c -Wall -std=c++0x
LDFLAGS=
COMMON_SOURCES=Alignment.cpp Parser.cpp ScoreMatrix.cpp IntronStorage.cpp Kernel.cpp
TARGET_SOURCES=main.cpp
TEST_SOURCES=test/t_parser.cpp test/tests.cpp test/t_matrix.cpp
COMMON_OBJECTS=$(COMMON_SOURCES:.cpp=.o)
TARGET_OBJECTS=$(TARGET_SOURCES:.cpp=.o)
TEST_OBJECTS=$(TEST_SOURCES:.cpp=.o)
EXECUTABLE=prosplign_parser
TEST_EXECUTABLE=test/t_prosplign_parser

.PHONY: test all target clean

all: target

target: $(EXECUTABLE)

test: $(TEST_EXECUTABLE)

# pull in dependency info for *existing* .o files
-include $(COMMON_OBJECTS:.o=.d)
-include $(TARGET_OBJECTS:.o=.d)
-include $(TEST_OBJECTS:.o=.d)

$(EXECUTABLE): $(COMMON_OBJECTS) $(TARGET_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

$(TEST_EXECUTABLE): $(COMMON_OBJECTS) $(TEST_OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@

.cpp.o: 
	$(CC) $(CFLAGS) $< -o $@
	$(CC) -MM $(CFLAGS) $< > $*.d

clean:
	rm -rf $(COMMON_OBJECTS) $(TEST_OBJECTS) $(TARGET_OBJECTS) $(EXECUTABLE) $(TEST_EXECUTABLE) *.d test/*.d
