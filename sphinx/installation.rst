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

**—————————— This section to be needed when pybind11 is used.**

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

**—————————— End of this section.**

You need to define the location of the dynamic libraries,
for ``ctypes`` to be able to import it. This is done through the ``pypw85.cfg``
file, which you must create and place in the following directory

- Windows 10/8/7/Vista: ``C:\Users\<username>\AppData\Roaming\pypw85``
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
name*. All these configure opions can be retrieved from the output of
``cmake --install .``. For example, on a Windows machine, where the output was::

  $ cmake --install .
  -- Install configuration: ""
  -- Installing: C:/opt/pw85/include/pw85
  -- Installing: C:/opt/pw85/include/pw85/pw85.h
  -- Installing: C:/opt/pw85/include/pw85/pw85_legacy.h
  -- Installing: C:/opt/pw85/lib/libpw85.dll.a
  -- Installing: C:/opt/pw85/lib/libpw85.dll
  -- Installing: C:/opt/pw85/lib/cmake/pw85/pw85-targets.cmake
  -- Installing: C:/opt/pw85/lib/cmake/pw85/pw85-targets-noconfig.cmake
  -- Installing: C:/opt/pw85/lib/cmake/pw85/pw85-config.cmake

the contents of ``pw85.ini`` is::

  [pw85]
  libpw85 = C:/opt/pw85/lib/libpw85.dll
  libpw85_legacy = C:/opt/pw85/lib/libpw85_legacy.dll
  datadir = C:/opt/pw85/share/pw85

To run the tests with Pytest_::

  $ python -m pytest tests/test_pw85.py

You can also test the “legacy” API. This requires the h5py_ module. To run the
tests, issue the command::

  $ python -m pytest tests/test_pw85_legacy.py

(beware, these tests take some time!).

.. _Breathe: https://breathe.readthedocs.io/
.. _CMake: https://cmake.org/
.. _Doxygen: https://www.doxygen.nl/
.. _Pytest: https://docs.pytest.org/
.. _Sphinx: https://www.sphinx-doc.org/
.. _h5py: https://www.h5py.org/

.. Local Variables:
.. fill-column: 80
.. End:
