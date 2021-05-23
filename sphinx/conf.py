# -*- coding: utf-8 -*-
import datetime
import pathlib

try:
    import pypw85

    project = pypw85.metadata["name"]
    author = pypw85.metadata["author"]
    copyright = pypw85.metadata["year"] + ", " + pypw85.metadata["author"]
    release = pypw85.metadata["version"]
    url = pypw85.metadata["url"]
except ImportError:
    from sphinx.util import logging

    logger = logging.getLogger(__name__)
    logger.warn("Could not import python module; some metadata is missing")
    project = "**MISSING PROJECT NAME**"
    copyright = "**MISSING COPYRIGHT**"
    author = "**MISSING AUTHOR**"
    release = "**MISSING VERSION**"
    url = "**MISSING URL**"

title = 'Documentation of the {} library'.format(project)
basename = project.lower()

with open('../README.md', mode='r', encoding='utf-8') as f:
    f.readline()
    f.readline()
    short_description = f.readline()

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.todo',
              'sphinx.ext.viewcode',
              'sphinx.ext.githubpages',
              'breathe']

numfig = True

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'

html_theme = 'sphinxdoc'
html_title = title
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

breathe_projects_source = {
    "pw85": (
        "../include/pw85",
        ["pw85.hpp"],
    )
}

breathe_doxygen_config_options = {"GENERATE_TODOLIST": "YES"}
