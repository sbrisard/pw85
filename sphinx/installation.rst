.. _installation:

************
Installation
************

.. contents:: Contents
   :local:

.. highlight:: none

Installing the C library
========================

The C library depends on

1. The `GLib <https://developer.gnome.org/glib/>`_ (for testing purposes)
2. The `GNU Scientific Library (GSL) <https://www.gnu.org/software/gsl/>`_ (for
   the implementation of the Brent algorithm)
3. The `HDF5 <https://portal.hdfgroup.org/>`_ (for testing purposes)

This is a CMake_ based project. The installation procedure is standard. First,
clone the repository. Then, ``cd`` into the root directory of the
pw85 project. Let
``pw85_INSTALL_PREFIX`` be the path to the directory
where pw85 should be installed::

  $ git clone https://github.com/sbrisard/pw85.git
  $ cd pw85
  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=pw85_INSTALL_PREFIX ..
  $ cmake --build . --config Release
  $ cmake --install . --config Release

.. note:: The ``--config`` option might not be available, depending on the
   selected generator.

At this point, pw85 should be installed. You can now
run the tests::

  $ ctest . -C Release

.. note:: Depending on the system, you might need to add
   ``pw85_INSTALL_PREFIX`` to your ``PATH`` environment
   variable.

To further test your installation, build the example in the :ref:`c-tutorial`.


Building the documentation
==========================

The documentation of pw85 requires Sphinx_. The C++ API
docs are built with Doxygen_ and the Breathe_ extension to Sphinx_.

To build the HTML version of the docs in the ``public`` subdirectory::

  $ cd docs
  $ sphinx-build -b html . ../public

To build the LaTeX version of the docs::

  $ cd docs
  $ make latex


Installing the Python bindings
==============================

To install the pypw85 module, ``cd`` into the
``python`` subdirectory and edit the ``setup.cfg`` file. Set the ``include_dir``
and ``library_dir`` to the appropriate paths. These should be::

  [pypw85]
  include_dir = ${CMAKE_INSTALL_PREFIX}/include
  library_dir = ${CMAKE_INSTLAL_PREFIX}/lib

Then, issue the following command::

  $ python setup.py install --user

or (if you intend to edit the project)::

  $ python setup.py develop --user

To run the tests with Pytest_::

  $ python -m pytest tests

.. todo:: Rewrite from there

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
``pypw85.cfg`` file should be::

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


.. _Breathe: https://breathe.readthedocs.io/
.. _CMake: https://cmake.org/
.. _Doxygen: https://www.doxygen.nl/
.. _Pytest: https://docs.pytest.org/
.. _Sphinx: https://www.sphinx-doc.org/
