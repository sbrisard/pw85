.. _installation:

************
Installation
************

.. contents:: Contents
   :local:

This chapter describes how to install the C library as well as the Python
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

Open the command prompt and ``cd`` into the ``src/build/`` directory. Issue the
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

Congratulations, ``PW85`` is now built and installed!

If you intend to use ``PW85`` from within Python only, then go to
:ref:`installation-of-the-python-bindings`.

If you also intend to link the library to e.g. C executables, you must inform
your system about the location of the library.

- On Windows platforms, you need to add the full path to ``pw85.dll`` to your
  ``PATH`` environment variable.
- On Linux or MacOS platforms, no further operation is required.

To test your installation, build the example in the :ref:`c-tutorial`.

.. _installation-of-the-python-bindings:

Installation of the Python bindings
===================================

The installation procedure is fairly standard and should be platform
independent. Installation in a virtual environment is not covered here, but is
possible with little alterations to the procedure below.

Open a terminal and ``cd`` into the ``wrappers/python-ctypes/``
directory. Issue the following command::

  $PYTHON_EXEC setup.py install

where ``$PYTHON_EXEC`` denotes your Python 3 executable (usually, ``python`` or
``python3``). Then, you need to define the location of the dynamic library, for
``ctypes`` to be able to import it. This is done through the ``pw85.ini`` file,
which you must create and place at the root of your home directory.

The contents of this file should be::

  [pw85]
  FullPath = full/path/to/the/pw85/dynamic/library

where the ``FullPath`` entry is the full path to the dynamic library
(``*.dll``, ``*.so`` or ``*.dylib``) *including its name*. It can be retrieved
from the output of ``cmake install``. For example, on a Windows machine, where
the output was::

  [100%] Built target pw85
  Install the project...
  -- Install configuration: "Release"
  -- Installing: C:/opt/pw85/lib/pw85-1.0/pw85-config.cmake
  -- Installing: C:/opt/pw85/lib/pw85-1.0/pw85-config-version.cmake
  -- Installing: C:/opt/pw85/lib/libpw85.dll.a
  -- Installing: C:/opt/pw85/bin/libpw85.dll
  -- Installing: C:/opt/pw85/lib/pw85-1.0/pw85-targets.cmake
  -- Installing: C:/opt/pw85/lib/pw85-1.0/pw85-targets-release.cmake
  -- Installing: C:/opt/pw85/include/pw85.h

the contents of ``pw85.ini`` is::

  [pw85]
  FullPath = C:/opt/pw85/bin/libpw85.dll

Provided the `pytest <https://pytest.org/>`_ module is installed on your
machine, you can run the tests as follows (from the ``wrappers/python-ctypes``
drectory)::

  $PYTHON_EXEC -m pytest tests

Beware! Tests take a while to run…
