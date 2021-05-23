import pathlib
import re

import pybind11
import setuptools

from setuptools.command.build_ext import build_ext


def read_metadata():
    metadata = {}
    filename = pathlib.Path.cwd() / ".." / "include" / "pw85" / "pw85.hpp"
    with open(filename, "r", encoding="utf8") as f:
        lines = f.readlines()
    prog = re.compile(
        r"constexpr\s*std::string_view\s*([a-z]*)\s*\{\s*\"([^\"]*)\"\s*}\s*;"
    )
    for line in lines:
        result = prog.match(line)
        if result is not None:
            metadata[result.group(1)] = result.group(2)
    return metadata


class my_build_ext(build_ext):
    def build_extensions(self):
        extra_compile_args = {
            "msvc": ["/std:c++17"],
            "unix": ["-std=c++17"],
        }

        for extension in self.extensions:
            extension.extra_compile_args += extra_compile_args[
                self.compiler.compiler_type
            ]
        build_ext.build_extensions(self)


if __name__ == "__main__":
    metadata = read_metadata()
    print(metadata)

    include_dirs = [
        pybind11.get_include(),
        pathlib.Path.cwd() / ".." / "include",
    ]

    pyfftwpp = setuptools.Extension(
        "pypw85",
        include_dirs=include_dirs,
        sources=["pypw85.cpp"],
        language="c++",
    )

    setuptools.setup(
        long_description_content_type="text/markdown",
        packages=setuptools.find_packages(),
        ext_modules=[pyfftwpp],
        cmdclass={"build_ext": my_build_ext},
        **metadata
    )
