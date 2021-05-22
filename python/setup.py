import pathlib

import pybind11
import setuptools

from setuptools.command.build_ext import build_ext


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
    include_dirs = [
        pybind11.get_include(),
        pathlib.Path.cwd() / ".." / "include",
        # TODO Path to boost should be a parameter
       r"C:\Users\sbrisard\miniconda3\Library\include"
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
    )
