.. _installation:

************
Installation
************

.. contents:: Contents
   :local:

This chapter describes how to install the C library as well as the
Python bindings. The first step is to clone the source of this
library. You can either get the latest release at
https://github.com/sbrisard/pw85/releases or clone the Git repository
for the development version::

  git clone https://github.com/sbrisard/pw85.git

.. highlight:: none

Building and installing the C library
=====================================

PW85 requires a POSIX system. On Windows platforms, it is recommended that you
install `MSYS2 <https://www.msys2.org/The>`_.

.. note:: Meson supports MSVC compilers. However, at the time of writing,
          MSVC-based installation has not been tested. Contributions are most
          welcome!

The C library depends on

1. The `GLib <https://developer.gnome.org/glib/>`_ (for testing purposes)
2. The `GNU Scientific Library (GSL) <https://www.gnu.org/software/gsl/>`_ (for
   the implementation of the Brent algorithm)
3. The `HDF5 <https://portal.hdfgroup.org/>`_ (for testing purposes)

The installation procedure also requires Python 3.

For installation, we use the `Meson build system
<https://mesonbuild.com/>`_. We assume that GLib, GSL and Meson are installed
on your system. Assuming that the project is built in the ``src/build/``
directory (no need to create it first), here is the whole installation
procedure (you must first ``cd`` into the root directory of the PW85 project)::

  cd pw85/src
  meson build
  cd build
  ninja install

A prefix can be specified in order for PW85 to be installed in a custom
directory, like so::

  cd pw85/src
  meson build --prefix=C:/opt/pw85
  cd build
  ninja install

Congratulations, ``PW85`` is now built and installed! You can then test the
installation (stay in the ``src/build`` directory)::

  meson test

If you intend to use ``PW85`` from within Python only, then go to
:ref:`installation-of-the-python-bindings`.

If you also intend to link the library to e.g. C executables, you must inform
your system about the location of the library.

- On Windows platforms, you need to add the full path to ``pw85.dll`` to your
  ``PATH`` environment variable.
- On Linux or MacOS platforms, no further operation is required.

To further test your installation, build the example in the :ref:`c-tutorial`.

.. _installation-of-the-python-bindings:

Installation of the Python bindings
===================================

The installation procedure is fairly standard and should be platform
independent. It requires a fairly recent version of `NumPy
<https://numpy.org/>`_. I successfully installed the Python bindings alongside
v1.16.4 of NumPy. Please report if you are successful with older versions.

Installation in a virtual environment is not covered here, but is possible with
little alterations to the procedure below.

Open a terminal and ``cd`` into the ``wrappers/python-ctypes/`` directory. Issue
the following command::

  $PYTHON_EXEC setup.py install

where ``$PYTHON_EXEC`` denotes your Python 3 executable (usually, ``python`` or
``python3``). Then, you need to define the location of the dynamic libraries,
for ``ctypes`` to be able to import it. This is done through the ``pypw85.cfg``
file, which you must create and place in the following directory

- Windows 10/8/7/Vista: ``C:\Users\<User Name>\AppData\Roaming\pypw85``
- Windows XP/2000: ``C:\Documents and Settings\<username>\Application
  Data\pypw85``
- Mac: ``/Users/<username>/Library/Application Support/pypw85``
- Linux: ``~/.pypw85``

(in all cases, the ``pypw85`` subdirectory must be created). The contents of the
``pw85.cfg`` file should be::

  [pw85]
  libpw85 = full/path/to/the/pw85/dynamic/library
  libpw85_legacy = full/path/to/the/pw85_legacy/dynamic/library
  datadir = full/path/to/the/pw85/data/directory

where the ``libpw85`` and ``libpw85_legacy`` entries are the full path to the
dynamic libraries (``*.dll``, ``*.so`` or ``*.dylib``) *including their
name*. All these configure opions can be retrieved from the output of ``ninja
install``. For example, on a Windows machine, where the output was::

  Installing libpw85.dll to C:/opt/pw85/bin
  Installing libpw85.dll.a to C:/opt/pw85/lib
  Installing libpw85.a to C:/opt/pw85/lib
  Installing libpw85_legacy.dll to C:/opt/pw85/bin
  Installing libpw85_legacy.dll.a to C:/opt/pw85/lib
  Installing libpw85_legacy.a to C:/opt/pw85/lib
  Installing pw85_ref_data.h5 to C:/opt/pw85/share/pw85
  Installing C:\path\to\pw85\src\pw85_legacy.h to C:/opt/pw85/include
  Installing C:\path\to\pw85\src\build\pw85.h to C:/opt/pw85/include

the contents of ``pw85.ini`` is::

  [pw85]
  libpw85 = C:/opt/pw85/bin/libpw85.dll
  libpw85_legacy = C:/opt/pw85/bin/libpw85_legacy.dll
  datadir = C:/opt/pw85/share/pw85

Provided the `pytest <https://pytest.org/>`_ module is installed on your
machine, you can run the tests as follows (from the ``wrappers/python-ctypes``
drectory)::

  $PYTHON_EXEC -m pytest tests/test_pw85.py

You can also test the “legacy” API. This requires the `h5py
<https://www.h5py.org/>`_ module. To run the tests, issue the command::

  $PYTHON_EXEC -m pytest tests/test_pw85_legacy.py

(beware, these tests take some time!).
