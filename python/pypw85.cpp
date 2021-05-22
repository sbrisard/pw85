#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pw85/pw85.hpp>

PYBIND11_MODULE(pypw85, m) {
  pybind11::dict metadata;
  metadata["author"] = pybind11::cast(pw85::metadata::author);
  metadata["description"] = pybind11::cast(pw85::metadata::description);
  metadata["author_email"] = pybind11::cast(pw85::metadata::author_email);
  metadata["license"] = pybind11::cast(pw85::metadata::license);
  metadata["name"] = pybind11::cast(pw85::metadata::name);
  metadata["url"] = pybind11::cast(pw85::metadata::url);
  metadata["version"] = pybind11::cast(pw85::metadata::version);
  metadata["year"] = pybind11::cast(pw85::metadata::year);
  m.attr("metadata") = metadata;

  m.def(
      "_cholesky_decomp",
      [](pybind11::array_t<double> a, pybind11::array_t<double> l) {
        pw85::_cholesky_decomp(a.data(), l.mutable_data());
      },
#include "docstrings/_cholesky_decomp.txt"
      ,
      pybind11::arg("a"), pybind11::arg("l"));
}
