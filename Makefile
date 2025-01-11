CC = gcc
CFLAGS = -g -Wall -O3
INCLUDES = -Isrc -Isrc/bgg -Isrc/sampling -Isrc/utils -Isrc -Isrc/cprf

SRC_DIR = src
BUILD_DIR = build
SRC_TEST = tests

OBJS_RAW = common matrix attribute random sampling circuit gen_circuit bgg cp cprf
OBJS_O = $(addsuffix .o,$(OBJS_RAW))
OBJS = $(addprefix $(BUILD_DIR)/,$(OBJS_O))

# list of executables binaries
EXEC_RAW = sampling circuit bgg cp_bit gen_circuit is_short cp kbitprf_circuit  kbitprf eval eval_circuit constrain constrain_circuit constrain_eval
EXEC = $(addprefix test_,$(EXEC_RAW))

# build cp as a library, including math lib
SPECIFY_LIBS = -L./build '-Wl,-rpath,./build'
LIBS_FLAGS = $(SPECIFY_LIBS) -lcp -lm


default: libcp.so

tests: $(EXEC)


# Building shared library
libcp.so: $(OBJS)
	$(CC) -shared -o $(BUILD_DIR)/$@ $^

# Brute-force searching
# Using -fPIC as the files will be included in shared lib
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -c -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/utils/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -c -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/bgg/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -c -o $@ $^
	
$(BUILD_DIR)/%.o: $(SRC_DIR)/cprf/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -c -o $@ $^

# -maes required for random.c using aes intrinsics
$(BUILD_DIR)/%.o: $(SRC_DIR)/sampling/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -maes -fPIC -c -o $@ $^


# Rule for building tests
test_%: $(SRC_TEST)/test_%.c
	$(CC) $(CFLAGS)  $(INCLUDES) -o $(BUILD_DIR)/$@ $^ $(LIBS_FLAGS)


# We clean object files and binaries
clean:
	rm -f $(BUILD_DIR)/*.o $(addprefix $(BUILD_DIR)/,$(EXEC))

# We clean the lib too
cleanall:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.so $(addprefix $(BUILD_DIR)/,$(EXEC))

.PHONY: tests clean cleanall
