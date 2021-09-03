# CHANGE HERE IF YOU INSTALLED CPLEX IN A NON-DEFAULT LOCATION
# (OR INSTALLED ANOTHER VERSION) :
# ---> CPLEX=/path/to/your/CPLEX/folder
CPLEX=/opt/ibm/ILOG/CPLEX_Studio1210

# project name (generate executable with this name)
TARGET   = biorseo
CC	   = g++
CFLAGS   = -Icppsrc/ -I/usr/local/include -I$(CPLEX)/concert/include -I$(CPLEX)/cplex/include -g -O3
CXXFLAGS = --std=c++17 -Wall -Wpedantic -Wextra -Wno-deprecated-copy -Wno-ignored-attributes
LINKER   = g++
LDFLAGS  = -L$(CPLEX)/concert/lib/x86-64_linux/static_pic/ -L$(CPLEX)/cplex/lib/x86-64_linux/static_pic/ -lboost_system -lboost_filesystem -lboost_program_options -lgomp -lconcert -lilocplex -lcplex -lpthread -ldl -lRNA -lm

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
	@echo -e "\033[00;32mLinking completed.\033[00m"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp $(INCLUDES)
	@mkdir -p $(OBJDIR)
	$(CC) -c $(CFLAGS) $(CXXFLAGS) $< -o $@
	@echo -e "\033[00;32mCompiled "$<".\033[00m"

doc: mainpdf supppdf
	@echo -e "\033[00;32mLaTeX documentation rendered.\033[00m"

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
	@echo -e "\033[00;32mCleanup completed.\033[00m"

.PHONY: remove
remove:
	@$(rm) $(BINDIR)/$(TARGET)
	@$(rm) doc/main_bioinformatics.pdf doc/supplementary_material.pdf
	@echo -e "\033[00;32mExecutable and docs removed!\033[00m"
