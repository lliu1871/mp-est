$(info $$PREFIX is [${PREFIX}])
# Report the PATH prefix

TARGET=mpest
$(info $$TARGET is [${TARGET}])
# Report the target

CC=gcc
$(info $$CC is [${CC}])
# Report the compiler used

# End of user configuration
############

CFLAGS=-DUNIX_VERSION -O3 -Wall -Wno-uninitialized
LDFLAGS= -lm
INCLUDES=src/mpest.o src/tool.o src/neighbour_joining.o
# Compile options

$(TARGET): $(INCLUDES)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(LDFLAGS)
# Compile the program

############

.PHONY: install
install: $(TARGET)
	cp $(TARGET) $(PREFIX)/bin/$(TARGET)
# Command to install by moving binary

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)
# Command to uninstall by removing binary

.PHONY: clean
clean:
	rm -f src/*.o $(TARGET)
# Command to remove all compiled files to make a clean install

############