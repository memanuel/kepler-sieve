# Modules for this project:
# None

# Tools and libraries used
# gcc version 10.3.0
# Boost 1.71
# (both gcc and Boost current version installed through apt)
# MariaDB Connector/C version 3.1.11 and Connector/C++ 1.0
# MariaDB server is on version 10.5.10

# Documentation for MariaDB Connector/C++:
# https://mariadb.com/docs/clients/connector-cpp/

# *****************************************************************************
# Compiler settings
CXX=g++
CFLAGS=-std=c++17

# Linker settings
LIBPATH=-L/usr/lib
LIBS=-lmariadbcpp
LDFLAGS=$(LIBPATH) $(LIBS)

# *****************************************************************************
# Targets
db_test: db_test.cpp db_utils.cpp
	$(CXX) db_test.cpp db_utils.cpp -o db_test.x $(CFLAGS) $(LDFLAGS)

mariadb-example: mariadb-example.cpp
	$(CXX) mariadb-example.cpp -o mariadb-example.x $(CFLAGS) $(LDFLAGS)

# *****************************************************************************
# Special targets
all: db_test

.PHONY: clean

clean:
	rm -f *.x
