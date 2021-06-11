This repository contains the plugin Affine for the Ibex library.

# Compilation of the plugin Affine

```
mkdir build
cd build
cmake ..
make
make check
```

If Ibex was installed in a directory unknown to CMake, you must set the
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
`find_package (ibex-affine)` in a CMake script.

To use Affine Form, first look at the file `examples/ex_affineform.cpp`. 
There is some simple examples to know how to use it.

To compile it, you can use the `makefile` and `pkg-config`. 
If you have specified some special path to install `ibex-affine` or `ibex-lib`, you should add the directory containing `ibex-affine.pc` and `ibex.pc`
to the PKG_CONFIG_PATH environment variable.

```
cd examples
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/my/special/path/to/ibex.pc/:/my/special/path/to/ibex-affine.pc/
make
./ex_affineform
```
