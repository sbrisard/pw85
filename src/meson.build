project('PW85', 'c', default_options: ['default_library=both'])

math = meson.get_compiler('c').find_library('m', required : false)
gsl = dependency('gsl')
glib = dependency('glib-2.0')
hdf5_hl = dependency('hdf5_hl')


# ------------
# CORE LIBRARY
# ------------

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


# --------------
# LEGACY LIBRARY
# --------------

install_headers('pw85_legacy.h')

pw85_legacy_lib = library('pw85_legacy',
			  ['pw85_legacy.c'],
			  link_with: pw85_lib,
			  dependencies: [math, gsl],
			  install: true)


# -----
# TESTS
# -----

src_dir = meson.current_source_dir()

pw85_ref_data_basename = 'pw85_ref_data.h5'
pw85_ref_data_path = join_paths(src_dir, '..', 'data', pw85_ref_data_basename)

python = import('python').find_installation('python3')

pw85_ref_data = custom_target(pw85_ref_data_basename,
			      output: pw85_ref_data_basename,
			      command: [python,
					join_paths(src_dir,
						   'download_ref_data.py'),
					'--pw85_ref_data_path',
					pw85_ref_data_path])

pw85_test = executable('test_pw85', 'test_pw85.c',
		       c_args: ['-DPW85_REF_DATA_PATH="'+pw85_ref_data_path+'"'],
		       link_with: pw85_lib,
		       dependencies: [glib, gsl, hdf5_hl])

test('GTest tests of PW85', pw85_test, depends: pw85_ref_data)
