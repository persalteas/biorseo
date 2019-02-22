include EditMe
ICONCERT=$(CPLEXDir)/concert/include
ICPLEX=$(CPLEXDir)/cplex/include
LCONCERT=$(CPLEXDir)/concert/lib/x86-64_linux/static_pic/
LCPLEX=$(CPLEXDir)/cplex/lib/x86-64_linux/static_pic/

# project name (generate executable with this name)
TARGET   = biominserter

CC	   = clang++
# compiling flags here
CFLAGS   = -Icppsrc/ -I/usr/local/include -I$(ICONCERT) -I$(ICPLEX) -I$(INUPACK) -I$(IEIGEN) -O3  
CXXFLAGS = --std=c++17 -Wall -Wpedantic -Wextra -Wno-ignored-attributes -Wno-unused-variable

LINKER   = clang++
# linking flags here
LDFLAGS   = -L$(LCONCERT) -L$(LCPLEX) -lboost_system -lboost_filesystem -lconcert -lilocplex -lcplex -lpthread -ldl -lnupackpfunc -lnupackutils

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
