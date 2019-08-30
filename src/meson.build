project('PW85', 'c', default_options: ['default_library=both'])

# Find dependencies
math = meson.get_compiler('c').find_library('m', required : false)
gsl = dependency('gsl')
glib = dependency('glib-2.0')

python = import('python').find_installation('python3')


pw85_header = custom_target('pw85.h',
			    input: ['preprocess_headers.py',
				    '../metadata/version.txt',
				    'pw85.h.in'],
			    output: 'pw85.h',
			    command: [python, '@INPUT@', '@OUTPUT@'],
			    install: true,
			    install_dir: 'include')

pw85_lib = library('pw85',
		   ['pw85.c', pw85_header],
		   dependencies: [math, gsl],
		   install: true)

pw85_test = executable('test_pw85', 'test_pw85.c',
		       link_with: pw85_lib,
		       dependencies: [glib, gsl])

install_headers('pw85_legacy.h')

pw85_legacy_lib = library('pw85_legacy',
			  ['pw85_legacy.c', 'pw85_legacy.h', pw85_header],
			  link_with: pw85_lib,
			  dependencies: [math, gsl],
			  install: true)

test('GTest tests of PW85', pw85_test)
