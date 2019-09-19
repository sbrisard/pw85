project('PW85_C_Tutorial', 'c')

pw85_include = include_directories(get_option('pw85_include'))

cc = meson.get_compiler('c')
pw85_lib = cc.find_library('pw85', dirs: get_option('pw85_lib'))

exe = executable('tutorial',
		 'tutorial.c',
		 include_directories: pw85_include,
		 dependencies: pw85_lib)
