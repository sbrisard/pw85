.. _installation:

************
Installation
************

This section describes how to install the C library as well as the Python
bindings. The first step is to clone the Git repository::

  git clone https://github.com/sbrisard/pw85.git

.. highlight:: none

Building and installing the C library
=====================================

The C library has no dependency per se. However, regardless of the platform,
the building process relies on `CMake <https://cmake.org/>`_, which you are
required to install. The instructions below use the command-line
exclusively. However, you can of course reach the same results from the
interfaces ``ccmake`` or ``cmake-gui``.

.. note:: The minimum required version of CMake is 3.0. Older versions of CMake
          have not been tested. Please report any success!

For all platforms, we will assume that the project is built in the
``src/build/`` directory (that you should first create)::

  cd pw85/src
  mkdir build
  cd build

Windows + Visual C++ platforms
------------------------------

Open the ``Visual C++ 2015 x64 Native Buld Tools Command Prompt`` (or
equivalent). This ensures that all native build tools will be correctly
discovered by CMake. ``cd`` into the ``src/build/`` directory. Issue the
following call to cmake::

  C:\path\to\pw85\src\build>C:\path\to\cmake\cmake.exe .. -G "NMake Makefiles" -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_INSTALL_PREFIX:PATH=C:/opt/pw85

(feel free to modify the install prefix). Then, build and install the project::

  C:\path\to\pw85\src\build>nmake
  C:\path\to\pw85\src\build>nmake install

Go to :ref:`finalize-your-installation`.

Windows + MinGW platforms
-------------------------

Open the command prompt and `cd`` into the ``src/build/`` directory. Issue the
following call to cmake::

  C:\path\to\cmake\bin\cmake.exe .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_INSTALL_PREFIX:PATH=C:/opt/pw85

(feel free to modify the install prefix). Then, build and install the project::

  C:\path\to\pw85\src\build>mingw32-make
  C:\path\to\pw85\src\build>mingw32-make install

Go to :ref:`finalize-your-installation`.

Linux platforms
---------------

Open a terminal and ``cd`` into the ``src/build/`` directory. Issue the
following call to cmake::

  cmake .. -G "Unix Makefiles" -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_INSTALL_PREFIX:PATH=~/local

(feel free to modify the install prefix). Then, build and install the project::

  make
  make install

Go to :ref:`finalize-your-installation`.

MacOS + HomeBrew platforms
--------------------------

.. _finalize-your-installation:

Finalize your installation
--------------------------

Congratulations, ``PW85`` is now built and installed! In order to make sure
that the dynamic library will be found in any circumstances (including from
within Python), you need to inform your system:

- on Windows platforms, you need to add the full path to ``pw85.dll`` to your ``PATH`` environment variable.

To test your installation, build the example in the :ref:`c-tutorial`.

Installation of the Python bindings
===================================
