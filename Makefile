# Compiler to be used
CC = gcc

# Directories
SRC_DIR = src
LIB_DIR = $(SRC_DIR)/lib
OBJ_DIR = obj

# Compiler flags - Added -I to specify include paths for header files
CFLAGS = -g -fopenmp -I$(SRC_DIR) -I$(LIB_DIR)

# Libraries to link against
LIBS = -lm -lgsl -lgslcblas

# The target executable
TARGET = RQA_in_RPS_systems

# Automatically find all .c files in the source directories
SRCS = $(wildcard $(SRC_DIR)/*.c) $(wildcard $(LIB_DIR)/*.c)

# Create a list of object files that will be placed in the OBJ_DIR
OBJS = $(patsubst %.c,$(OBJ_DIR)/%.o,$(notdir $(SRCS)))

# Default target to build and run
run: $(TARGET)
	@./$(TARGET)

# Rule to link object files from the obj/ directory to create the executable
$(TARGET): $(OBJS)
	@echo "==> Linking to create $(TARGET)..."
	@$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# This is a pattern rule that tells make how to build any .o file in the obj/ directory.
# It checks for the corresponding .c file in both src/ and src/lib/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	@echo "==> Compiling $<..."
	@$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(LIB_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	@echo "==> Compiling $<..."
	@$(CC) $(CFLAGS) -c $< -o $@

# Clean up the build (remove the obj directory and the executable)
clean:
	@echo "==> Cleaning up build files..."
	@rm -rf $(OBJ_DIR) $(TARGET)

# Phony targets to avoid conflicts with files of the same name
.PHONY: run clean