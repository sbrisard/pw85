.. _installation

************
Installation
************

This section describes how to install the C library as well as the Python
bindings.


Building and installing the C library
=====================================

The C library has no dependency per se. However, regardless of the platform,
the building process relies on `CMake <https://cmake.org/>`_, which you are
required to install. Depending on the platforms, the instructions may use the
command line tool ``ccmake`` or the GUI ``cmake-gui``.

For all platforms, we will assume that the project will be built in the
``src/build`` directory (that you should first create).

Linux platforms
---------------


Windows + Visual C++ platforms
------------------------------

Open the ``Visual C++ 2015 x64 Native Buld Tools Command Prompt`` (or
equivalent), and from there, launch ``cmake-gui.exe`` (this ensures that all
native build tools will be correctly discovered by CMake)::

  C:\path\to\cmake\bin\cmake-gui.exe

In the main window, specify “where is the source code” and “where to build the
binaries” (see below)

.. figure:: installation/windows-where.*
   :align: center

Click ``Configure``, and specify the “NMake Makefiles” generator, using default
native compilers (see below).

.. figure:: installation/windows-specify_nmake.*
   :align: center

Click ``Finish``. CMake should echo the “Configuring done” message (see
below). Set the variable ``CMAKE_BUILD_TYPE`` to ``Release`` (**this is
important!**); check the ``CMAKE_INSTALL_PREFIX``.

.. figure:: installation/windows-finish.*
   :align: center

Modify any options and click ``Configure`` again when you are done. Then click
``Generate``. Close GUI, and go back to the
``Visual C++ 2015 x64 Native Build Tools Command Prompt``. ``cd`` into the ``build/`` directory, and issue the following command::

  nmake install

This is what your console should look like

.. figure:: installation/windows-run_nmake.*
   :align: center

Congratulations, ``PW85`` is now properly built!

Windows + MinGW platforms
-------------------------


MacOS + HomeBrew platforms
--------------------------

Installation of the Python bindings
===================================
