This repository contains the plugin Affine for the Ibex library.

# Compilation of the plugin

```
mkdir build
cd build
cmake ..
make
make check
```

If Ibex was install in a directory unknown to CMake, you must set the
environment variable `IBEX_DIR` to indicate this directory, before calling
CMake. For example,

```
IBEX_DIR=/my/special/path cmake ..
```

# Installation

```
sudo make install
```

To install to a specific directory, you need to specify it when calling CMake.
For example,
```
  IBEX_DIR=/my/special/path cmake -DCMAKE_INSTALL_PREFIX=/another/path ..
```

# Using the plugin

To use the plugin, you can either use `pkg-config` or use
`find_package (Ibex-affine)` in a CMake script.
