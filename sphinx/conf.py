# -*- coding: utf-8 -*-
import datetime
import inspect

import pypw85

project = 'PW85'
title = 'Documentation of the {} library'.format(project)
author = 'Sébastien Brisard'
copyright = '2018–{}, {}'.format(datetime.date.today().year, author)
description = inspect.getdoc(pypw85).split('\n')[0]
version = pypw85.__version__
release = version
url = 'https://github.com/sbrisard/pw85.git'
basename = 'pw85'

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.todo',
              'sphinx.ext.viewcode',
              'sphinx.ext.githubpages',]

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'

html_theme = 'alabaster'
htmlhelp_basename = basename+'doc'

latex_engine = 'xelatex'
latex_elements = {'fontpkg': r'''
\setmainfont{Liberation Serif}
\setsansfont{Liberation Sans}
\setmonofont[Scale=MatchLowercase]{DejaVu Sans Mono}'''}
latex_documents = [(master_doc, basename+'.tex', title, author, 'manual'),]

man_pages = [(master_doc, basename, title, [author], 1)]

texinfo_documents = [(master_doc, basename, title, author, project, description,
                      'Miscellaneous'),]

epub_basename = basename
epub_title = title
epub_author = author
epub_publisher = author
epub_copyright = copyright
epub_identifier = url
epub_exclude_files = ['search.html']

todo_include_todos = True

autodoc_default_options = {'member-order': 'groupwise',
                           'undoc-members': True,}

# Create `abstract.rst` file
with open('abstract.rst', mode='w', encoding='utf-8') as f:
    f.write(inspect.getdoc(pypw85).replace('This module', project))
