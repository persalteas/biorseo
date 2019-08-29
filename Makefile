include EditMe
ICONCERT=/opt/ibm/ILOG/CPLEX_Studio_Community128/concert/include
ICPLEX=/opt/ibm/ILOG/CPLEX_Studio_Community128/cplex/include
INUPACK=/usr/local/include/nupack
IEIGEN=/usr/local/include/eigen3
LCONCERT=/opt/ibm/ILOG/CPLEX_Studio_Community128/concert/lib/x86-64_linux/static_pic/
LCPLEX=/opt/ibm/ILOG/CPLEX_Studio_Community128/cplex/lib/x86-64_linux/static_pic/

# project name (generate executable with this name)
TARGET   = biorseo

CC	   = clang++
# compiling flags here
CFLAGS   = -Icppsrc/ -I/usr/local/include -I$(ICONCERT) -I$(ICPLEX) -I$(INUPACK) -I$(IEIGEN) -O3
CXXFLAGS = --std=c++17 -Wall -Wpedantic -Wextra -Wno-ignored-attributes -Wno-unused-variable

LINKER   = clang++
# linking flags here
LDFLAGS   = -L$(LCONCERT) -L$(LCPLEX) -lboost_system -lboost_filesystem -lboost_program_options -lconcert -lilocplex -lcplex -lpthread -ldl -lnupackpfunc -lnupackutils

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

doc: mainpdf supppdf
	@echo "\033[00;32mLaTeX documentation rendered.\033[00m"

mainpdf: doc/main_bioinformatics.tex doc/references.bib doc/bioinfo.cls doc/natbib.bst
	cd doc; pdflatex -synctex=1 -interaction=nonstopmode -file-line-error main_bioinformatics
	cd doc; bibtex main_bioinformatics
	cd doc; pdflatex -synctex=1 -interaction=nonstopmode -file-line-error main_bioinformatics
	cd doc; pdflatex -synctex=1 -interaction=nonstopmode -file-line-error main_bioinformatics

supppdf: doc/supplementary_material.tex
	cd doc; pdflatex -synctex=1 -interaction=nonstopmode -file-line-error supplementary_material

.PHONY: all
all: $(BINDIR)/$(TARGET) doc

.PHONY: re
re: remove clean all

.PHONY: clean
clean:
	$(rm) $(OBJECTS)
	$(rm) doc/supplementary_material.bbl doc/supplementary_material.blg doc/supplementary_material.synctex.gz doc/supplementary_material.log doc/supplementary_material.aux
	$(rm) doc/main_bioinformatics.bbl doc/main_bioinformatics.blg doc/main_bioinformatics.synctex.gz doc/main_bioinformatics.log doc/main_bioinformatics.aux doc/OUP_First_SBk_Bot_8401-eps-converted-to.pdf
	@echo "\033[00;32mCleanup completed.\033[00m"

.PHONY: remove
remove:
	@$(rm) $(BINDIR)/$(TARGET)
	@$(rm) doc/main_bioinformatics.pdf doc/supplementary_material.pdf
	@echo "\033[00;32mExecutable and docs removed!\033[00m"
