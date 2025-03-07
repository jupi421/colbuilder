# Makefile for building the contact map tools for Colbuilder
# This Makefile should be placed in the topology directory

# Determine directory paths automatically
HOME ?= $(shell echo $$HOME)/bin
TOPO_DIR=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
FORCE_FIELD_DIR ?= $(TOPO_DIR)/../../topology
CONTACT_MAP_DIR=$(TOPO_DIR)
CONTACT_MAP=$(CONTACT_MAP_DIR)/contact_map
BIN_TARGET=$(HOME)/contact_map

# Compiler settings
CC?=clang -Wall
CFLAGS=-march=native -O3
LIBS=-lm

# Source files and dependencies
DEPS=$(CONTACT_MAP_DIR)/contact_map.h
SRC=$(CONTACT_MAP_DIR)/chemical_map.c $(CONTACT_MAP_DIR)/pdb_map.c $(CONTACT_MAP_DIR)/protein_map.c $(CONTACT_MAP_DIR)/contact_map.c
OBJ=$(SRC:.c=.o)

# Print all important paths before doing anything else
.PHONY: all print_paths check_contact_map rebuild_contact_map clean install

all: print_paths check_contact_map

print_paths:
	@echo "==== Path Configuration ===="
	@echo "HOME:             $(HOME)"
	@echo "TOPO_DIR:         $(TOPO_DIR)"
	@echo "FORCE_FIELD_DIR:  $(FORCE_FIELD_DIR)"
	@echo "CONTACT_MAP_DIR:  $(CONTACT_MAP_DIR)"
	@echo "CONTACT_MAP:      $(CONTACT_MAP)"
	@echo "BIN_TARGET:       $(BIN_TARGET)"
	@echo "DEPS:             $(DEPS)"
	@echo "SRC:              $(SRC)"
	@echo "OBJ:              $(OBJ)"
	@echo "============================"

# Check if contact map tool exists and build if needed
check_contact_map:
	@if [ ! -f "$(CONTACT_MAP)" ]; then \
		echo "Contact map tool not found, building..."; \
		$(MAKE) rebuild_contact_map; \
	else \
		echo "Contact map tool found at $(CONTACT_MAP)"; \
		if [ ! -f "./contact_map" ] && [ -d "." ]; then \
			ln -sf $(CONTACT_MAP) ./contact_map; \
			echo "Created symlink in ./"; \
		fi; \
	fi


# Force a rebuild of the contact map tool
rebuild_contact_map: clean $(CONTACT_MAP)
	@echo "Contact map tool built successfully at $(CONTACT_MAP)"

# Compile the contact map tool
$(CONTACT_MAP): $(OBJ)
	$(CC) -o $@ $^ $(LIBS)

# Compile object files
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Clean up compiled files
clean:
	rm -f $(CONTACT_MAP_DIR)/*.o
	rm -f $(CONTACT_MAP)
	rm -f $(BIN_TARGET)

# Install the contact map tool to the user's bin directory
install: check_contact_map
	@mkdir -p $(dir $(BIN_TARGET))
	@if [ -L "$(BIN_TARGET)" ]; then \
		echo "Removing old symbolic link..."; \
		rm $(BIN_TARGET); \
	fi
	@echo "Creating symbolic link to $(BIN_TARGET)..."
	ln -sf $(CONTACT_MAP) $(BIN_TARGET)
	@echo "Installation complete. The contact map tool is now available at: $(BIN_TARGET)"

# Helper to check dependencies
check_dependencies:
	@echo "Checking dependencies for Martini topology generation..."
	@echo "Contact map tool: $(if $(wildcard $(CONTACT_MAP)),Found,Missing - will be built)"
	@echo "create_goVirt.py: $(if $(wildcard $(TOPO_DIR)/create_goVirt.py),Found,Missing)"
	@command -v martinize2 >/dev/null 2>&1 && echo "Martinize2: Found" || echo "Martinize2: Missing - install with conda"
	@echo "ITP module: $(if $(wildcard $(TOPO_DIR)/../colbuilder/core/topology/itp.py),Found,Missing)"
