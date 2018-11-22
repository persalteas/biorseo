ICONCERT=/opt/ibm/ILOG/CPLEX_Studio_Community128/concert/include
ICPLEX=/opt/ibm/ILOG/CPLEX_Studio_Community128/cplex/include
LCONCERT=/opt/ibm/ILOG/CPLEX_Studio_Community128/concert/lib/x86-64_linux/static_pic/
LCPLEX=/opt/ibm/ILOG/CPLEX_Studio_Community128/cplex/lib/x86-64_linux/static_pic/

# project name (generate executable with this name)
TARGET   = biominserter

CC	   = g++
# compiling flags here
CFLAGS   = -Icppsrc/ -I$(ICONCERT) -I$(ICPLEX) -O3
CXXFLAGS = -std=c++17 -Wall -Wpedantic -Wextra -Wno-ignored-attributes

LINKER   = g++
# linking flags here
LDFLAGS   = -lconcert -lilocplex -lcplex -lm -lpthread -ldl -lboost_system -lboost_filesystem -lboost_program_options -L$(LCONCERT) -L$(LCPLEX)

# change these to proper directories where each file should be
SRCDIR   = cppsrc
OBJDIR   = obj
BINDIR   = bin

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm	   = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(LINKER) $(OBJECTS) $(LDFLAGS) -o $@
	@echo "\033[00;32mLinking completed.\033[00m"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(INCLUDES)
	@mkdir -p $(OBJDIR)
	$(CC) -c $(CFLAGS) $(CXXFLAGS) $< -o $@
	@echo "\033[00;32mCompiled "$<".\033[00m"

.PHONY: clean
clean:
	$(rm) $(OBJECTS)
	@echo "\033[00;32mCleanup completed.\033[00m"

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "\033[00;32mExecutable removed!\033[00m"