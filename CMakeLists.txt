# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.6 )
project( Cappucino )
# Locate Houdini's libraries and header files.
# Registers an imported library target named 'Houdini'.
find_package( Houdini REQUIRED )

# Add a library with two source files.
set( library_name Cappucino )
add_library( ${library_name} SHARED
    Main.C
    vdbConvolve.h
	vdbConvolve.C
	vdbCpt.h
	vdbCpt.C
	vdbWave.h
	vdbWave.C
	vdbWaveKernel.h
	vdbWaveKernel.C
	vdbDivergence.h
	vdbDivergence.C
	vdbReact.h
	vdbReact.C
)
# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini ${_houdini_root}/custom/houdini/dsolib/openvdb_sesi.lib ${_houdini_root}/custom/houdini/dsolib/half.lib)
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )
