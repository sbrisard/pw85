project('PW85', 'c', default_options: ['default_library=both'])

gsl = dependency('gsl')

prog_python = import('python').find_installation('python3')

pw85_h = custom_target('pw85.h',
		       input: ['preprocess_headers.py',
			       '../metadata/version.txt',
			       'pw85.h.in'],
		       output: 'pw85.h',
		       command: [prog_python, '@INPUT@', '@OUTPUT@'],
		       install: true,
		       install_dir: 'include')

pw85_lib = library('pw85',
		   ['pw85.c', pw85_h],
		   dependencies: gsl,
		   install: true)
