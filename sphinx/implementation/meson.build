project('f_accuracy', 'c')
glib = dependency('glib-2.0')
hdf5_hl = dependency('hdf5_hl')

pw85_include = include_directories(get_option('pw85_include'))

cc = meson.get_compiler('c')
pw85_lib = cc.find_library('pw85', dirs: get_option('pw85_lib'))
pw85_legacy_lib = cc.find_library('pw85_legacy', dirs: get_option('pw85_lib'))

pw85_ref_data_path = join_paths(get_option('pw85_data'), 'pw85_ref_data.h5')

exe = executable('f_accuracy',
		 'f_accuracy.c',
		 c_args: ['-DPW85_REF_DATA_PATH="'+pw85_ref_data_path+'"',
			  '-DHISTOGRAM_PATH="'+join_paths(meson.current_source_dir(), 'histograms.csv')+'"'],
		 include_directories: pw85_include,
		 dependencies: [glib, hdf5_hl, pw85_lib, pw85_legacy_lib])
