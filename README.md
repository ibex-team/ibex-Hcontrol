This repository contains a template to create plugins for the Ibex library with
CMake.

# Customization

To customize this template, replace all the occurences of `myplugin`,
`MYPLUGIN` and `MyPlugin` with the actual name of the plugin, and set all the
corresponding variables to their correct value.

If you are using the Git repository, you can use the following command to find
them all:
```
git grep -i myplugin
```

For example, in [CMakeLists.txt](./CMakeLists.txt), the two following lines
```
project (IBEX_MYPLUGIN VERSION 0.42 LANGUAGES CXX)
set (IBEX_MYPLUGIN_DESCRIPTION "Short description of the plugin")
```
could become
```
project (IBEX_AFFINE VERSION 1.0 LANGUAGES CXX)
set (IBEX_AFFINE_DESCRIPTION "Affine arithmetic for the Ibex library")
```


In [src/CMakeLists.txt](./src/CMakeLists.txt), there is an example of the
compilation of a library depending on the Ibex library.
In [src/bin/CMakeLists.txt](./src/bin/CMakeLists.txt), there is two examples of
the compilation of a binary.
One can remove the lines corresponding to examples that are not needed and can
customize the other ones (mostly by appending new source files and changing the
name of the targets).


The [CMakeLists.txt](./CMakeLists.txt) contains code to generate a pkg-config
file and CMake config file. If your plugin is installing a library, you must:
  - rename the two files `ibex-myplugin.pc.in` and
  `ibex-myplugin-config.cmake.in`;
  - change the paths accordingly in [CMakeLists.txt](./CMakeLists.txt);
  - check that all the variables and fields have the correct values.

If your plugin is not installing a library, you can:
  - remove the two files mentioned above;
  - remove in [CMakeLists.txt](./CMakeLists.txt) the code generating them.

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

Installation
------------

```
sudo make install
```

To install to a specific directory, you need to specify it when calling CMake.
For example,
```
  IBEX_DIR=/my/special/path cmake -DCMAKE_INSTALL_PREFIX=/another/path ..
```

# Dependencies

If your plugin requires Ibex to have been compiled with a specific library for
the interval arithmetic or the linear programming arithmetic, you can use the
variables `IBEX_INTERVAL_LIB_NAME` and `IBEX_LP_LIB_NAME`.

For example, if your plugin requires Ibex to have been compiled with a linear
programming library, you add in [CMakeLists.txt](./CMakeLists.txt):
```
# check if Ibex was compiled with an linear programming library
if (NOT IBEX_LP_LIB_NAME)
  message (FATAL_ERROR "A Linear Programming library is needed for this plugin")
endif ()
```
