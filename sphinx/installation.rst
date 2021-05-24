.. _installation:

************
Installation
************

.. contents:: Contents
   :local:

First of all, clone the repository

.. code-block:: none

  $ git clone https://github.com/sbrisard/pw85


Installing the C++ library
==========================

pw85 is a header-only library: there is no installation procedure *per se* and
you can drop the header wherever you like (as long as it is located in a
``pw85`` subdirectory). To use pw85 in a C++ project, you must include the
header


.. code-block:: cpp

   #include <pw85/pw85.hpp>

and inform the compiler of its location.

.. note:: pw85 depends on `Boost::Math <https://www.boost.org/doc/libs/1_75_0/libs/math/>`_
	  (for the implementation of the Brent algorithm). You must pass the
	  relevant options to the compiler. Typically, these would be ``-I``
	  options. The C++ tutorials provides a :ref:`CMake example <sec20210523205251>`.

To run the tests or build the documentation properly, you need to first build
the python bindings (see :ref:`below <sec20210523203528>`).


To further test your installation, build the example in the :ref:`sec20210523205251`.


.. _sec20210523203528:

Installing the Python bindings
==============================

The Python bindings are built with pybind11_, which must be installed.

To install the pw85 module, ``cd`` into the ``python`` subdirectory and run
the ``setup.py`` script as follows.

First, build the extension::

  $ python setup.py build_ext -Ipath/to/boost/math

When the extension is built, installation is down as usual::

  $ python setup.py install --user

or (if you intend to edit the project)::

  $ python setup.py develop --user

To run the tests with Pytest_::

  $ python -m pytest tests

(beware, these tests take some time!).

.. note:: Upon first execution, the test script will attempt to retrieve some
          precomputed reference data. In case of failure (e.g. if you sit behind
          a firewall), this reference file can be downloaded manually at this
          address: https://zenodo.org/record/3323683/files/pw85_ref_data-20190712.h5

	  The file should be placed in the ``data/`` subdirectory, at the root
	  of the project, and should be renamed ``pw85_ref_data.h5``::

              ├───data
              │   └───pw85_ref_data.h5
              ├───docs
              ├───include
              │   └───pw85
              ├───joss
              ├───legacy
              ├───python
              │   ├───docstrings
              │   └───tests
              └───sphinx
                  ├───cpp_tutorial
                  ├───implementation
                  │   └───f_accuracy
                  └───py_tutorial


Building the documentation
==========================

.. note:: For the documentation to build properly, the python module
          must be installed, as it is imported to retrieve the project
          metadata.

The documentation of pw85 requires Sphinx_. The C++ API docs are built with
Doxygen_ and the Breathe_ extension to Sphinx_.

To build the HTML version of the docs in the ``docs`` subdirectory::

  $ cd docs
  $ sphinx-build -b html . ../docs

To build the LaTeX version of the docs::

  $ cd docs
  $ make latex


.. _Breathe: https://breathe.readthedocs.io/
.. _CMake: https://cmake.org/
.. _Doxygen: https://www.doxygen.nl/
.. _pybind11: https://pybind11.readthedocs.io/
.. _Pytest: https://docs.pytest.org/
.. _Sphinx: https://www.sphinx-doc.org/
.. _h5py: https://www.h5py.org/

.. Local Variables:
.. fill-column: 80
.. End:
