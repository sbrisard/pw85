#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pw85/pw85.hpp>

using DoubleArray = pybind11::array_t<double>;

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

  m.attr("lambda_atol") = pw85::lambda_atol;

  m.def(
      "_cholesky_decomp",
      [](DoubleArray a, DoubleArray l) {
        pw85::_cholesky_decomp(a.data(), l.mutable_data());
      },
#include "docstrings/_cholesky_decomp.txt"
      , pybind11::arg("a"), pybind11::arg("l"));

  m.def(
      "_cholesky_solve",
      [](DoubleArray l, DoubleArray b, DoubleArray x) {
        pw85::_cholesky_solve(l.data(), b.data(), x.mutable_data());
      },
#include "docstrings/_cholesky_solve.txt"
      , pybind11::arg("l"), pybind11::arg("b"), pybind11::arg("x"));

  m.def(
      "spheroid",
      [](double a, double c, DoubleArray n, DoubleArray q) {
        pw85::spheroid(a, c, n.data(), q.mutable_data());
      },
#include "docstrings/spheroid.txt"
      , pybind11::arg("a"), pybind11::arg("c"), pybind11::arg("n"),
      pybind11::arg("q"));

  m.def(
      "f_neg",
      [](double lambda, DoubleArray r12, DoubleArray q1, DoubleArray q2) {
        return pw85::f_neg(lambda, r12.data(), q1.data(), q2.data());
      },
#include "docstrings/f_neg.txt"
      , pybind11::arg("lambda"), pybind11::arg("r12"), pybind11::arg("q1"),
      pybind11::arg("q2"));

  m.def(
      "contact_function",
      [](DoubleArray r12, DoubleArray q1, DoubleArray q2, DoubleArray out) {
        return pw85::contact_function(r12.data(), q1.data(), q2.data(),
                                      out.mutable_data());
      },
#include "docstrings/contact_function.txt"
      , pybind11::arg("r12"), pybind11::arg("q1"), pybind11::arg("q2"),
      pybind11::arg("out"));
}
