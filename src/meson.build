project('PW85', 'c', default_options: ['default_library=both'])

math = meson.get_compiler('c').find_library('m', required : false)
gsl = dependency('gsl')
glib = dependency('glib-2.0')

# Generate pw85.h header file (insert correct value for PW85_VERSION)
pw85_version = run_command('more', # `more` works on both windows and linux
			   join_paths(meson.source_root(),
				      '..',
				      'metadata',
				      'version.txt')).stdout().strip()

conf_data = configuration_data()
conf_data.set('project_version', pw85_version)
configure_file(input: 'pw85.h.in',
	       output: 'pw85.h',
	       configuration: conf_data,
	       install_dir: 'include')

pw85_lib = library('pw85', ['pw85.c'],
		   dependencies: [math, gsl],
		   install: true)

hdf5_hl = dependency('hdf5_hl')

pw85_test = executable('test_pw85', 'test_pw85.c',
		       link_with: pw85_lib,
		       dependencies: [glib, gsl, hdf5_hl])



test('GTest tests of PW85', pw85_test)

install_headers('pw85_legacy.h')

pw85_legacy_lib = library('pw85_legacy',
			  ['pw85_legacy.c'],
			  link_with: pw85_lib,
			  dependencies: [math, gsl],
			  install: true)


executable('hdf5_example',
	   'hdf5_example.c',
	   link_with: pw85_lib,
	   dependencies: [glib, hdf5_hl],
	   install: false)
