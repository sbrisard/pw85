import os
import distutils.log
import pathlib

from setuptools import setup
from setuptools.command.build_py import build_py

from string import Template

def get_metadata(basename, dir_):
    with open(dir_/(basename+'.txt'), mode='r', encoding='utf-8') as f:
        return f.read().strip()


def substitute(path, **kwds):
    with open(path, mode='r', encoding='utf-8') as infile, open(path+'.tmp', mode='w', encoding='utf-8') as outfile:
        for line in infile.readlines():
            outfile.write(Template(line).substitute(**kwds))
    os.replace(path+'.tmp', path)


class my_build_py(build_py):
    def run(self):
        super().run()
        metadata = self.distribution.metadata
        path = self.get_module_outfile(self.build_lib, ['pypw85'], '__init__')
        substitute(path,
                   author=metadata.get_author(),
                   version=metadata.get_version(),
                   long_description=metadata.get_long_description())


if __name__ == '__main__':
    name = 'pypw85'
    metadata_dir = pathlib.Path('.')/'..'/'..'/'metadata'
    version = get_metadata('version', metadata_dir)
    author = get_metadata('author', metadata_dir)
    url = get_metadata('url', metadata_dir)

    with open(pathlib.Path('../../README.md'), mode='r', encoding='utf-8') as f:
        lines = f.readlines()
        description = lines[2].strip()
        long_description = ''.join(lines[2:])

    setup(cmdclass={'build_py': my_build_py},
          packages=[name],
          name=name,
          version=version,
          author=author,
          author_email='',
          description=description,
          long_description=long_description,
          long_description_content_type='text/markdown',
          url=url,
          license='BSD-3',)
