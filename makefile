# Project Name (executable)
PROJECT = hypocoverage
# Compiler
CC = g++

# Run Options       
COMMANDLINE_OPTIONS = /dev/ttyS0

# Compiler options during compilation
COMPILE_OPTIONS =  -std=c++17 -pthread -O3 -w  -Ilibs -Iinclude

#Header include directories
HEADERS = 
#Libraries for linking
LIBS = -lpthread

# Dependency options
DEPENDENCY_OPTIONS =  -MM -std=c++17 -Ilibs -Iinclude






#-- Do not edit below this line --

# Subdirs to search for additional source files
SUBDIRS := $(shell ls src | grep "\/" )
DIRS := ./ $(SUBDIRS)
SOURCE_FILES := $(foreach d, $(DIRS), $(wildcard src/$(d)*.cpp) )

# Create an object file of every cpp file
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCE_FILES))

# Dependencies
DEPENDENCIES = $(patsubst %.cpp, %.d, $(SOURCE_FILES))

# Create .d files
%.d: %.cpp
	$(CC) $(DEPENDENCY_OPTIONS) $< -MT "$*.o $*.d" -MF $*.d

# Make $(PROJECT) the default target
all: $(DEPENDENCIES) $(PROJECT)

$(PROJECT): $(OBJECTS)
	$(CC) -o $(PROJECT) $(OBJECTS) $(LIBS)

# Include dependencies (if there are any)
ifneq "$(strip $(DEPENDENCIES))" ""
  include $(DEPENDENCIES)
endif

# Compile every cpp file to an object
%.o: %.cpp
	$(CC) -c $(COMPILE_OPTIONS) -o $@ $< $(HEADERS)

# Build & Run Project
run: $(PROJECT)
	./$(PROJECT) 

# Clean & Debug
.PHONY: makefile-debug
makefile-debug:

.PHONY: clean
clean:
	rm -f $(PROJECT) $(OBJECTS)

.PHONY: depclean
depclean:
	rm -f $(DEPENDENCIES)

clean-all: clean depclean
