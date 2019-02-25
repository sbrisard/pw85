# -*- coding: utf-8 -*-
import datetime
import pathlib


metadata_dir = pathlib.Path('../metadata')

def read_metadata(basename):
    with open(metadata_dir/(basename+'.txt'), 'r', encoding='utf-8') as f:
        return f.read().strip()

basenames = ['name', 'version', 'author', 'url', 'short_description']
project, version, author, url, short_description = [read_metadata(basename)
                                                    for basename in basenames]
release = version
copyright = '{}–{}, {}'.format(read_metadata('copyright'),
                               datetime.date.today().year,
                               author)

title = 'Documentation of the {} library'.format(project)
basename = project.lower()

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

texinfo_documents = [(master_doc, basename, title, author, project,
                      short_description, 'Miscellaneous'),]

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
