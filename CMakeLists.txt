cmake_minimum_required(VERSION 3.16)

project(minibwa VERSION 0.3 LANGUAGES C)

include(CheckCCompilerFlag)
include(GNUInstallDirs)

option(MINIBWA_ENABLE_GPL "Build GPL-licensed low-memory BWT construction sources" ON)
option(MINIBWA_ENABLE_MIMALLOC "Build with bundled mimalloc override" ON)
option(MINIBWA_ENABLE_OPENMP "Use OpenMP in bundled libsais when available" ON)
option(MINIBWA_BUILD_EXAMPLES "Build API examples" ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

set(MINIBWA_LIB_SOURCES
	src/kommon.c
	src/kalloc.c
	src/bwt.c
	src/l2bit.c
	src/options.c
	src/seed.c
	src/map-algo.c
	src/lchain.c
	src/align.c
	src/pe.c
	src/cs.c
	src/format.c
	src/ksw2_extz2_sse.c
	src/ksw2_extd2_sse.c
	src/ksw2_ll_sse.c
)

set(MINIBWA_CLI_SOURCES
	src/kthread.c
	third_party/libsais/libsais.c
	third_party/libsais/libsais64.c
	src/index.c
	src/bseq.c
	src/map-main.c
	src/fastmap.c
)

add_library(minibwa STATIC ${MINIBWA_LIB_SOURCES})
target_include_directories(minibwa
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_compile_options(minibwa PRIVATE -Wall)

add_executable(minibwa_cli src/main.c ${MINIBWA_CLI_SOURCES})
set_target_properties(minibwa_cli PROPERTIES OUTPUT_NAME minibwa)
target_include_directories(minibwa_cli PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/libsais
)
target_link_libraries(minibwa_cli PRIVATE minibwa Threads::Threads ZLIB::ZLIB m)
target_compile_options(minibwa_cli PRIVATE -Wall)

if(MINIBWA_ENABLE_GPL)
	target_sources(minibwa_cli PRIVATE
		third_party/bwtgen/QSufSort.c
		third_party/bwtgen/bwtgen.c
	)
	target_include_directories(minibwa_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/bwtgen)
	target_compile_definitions(minibwa PRIVATE USE_GPL)
	target_compile_definitions(minibwa_cli PRIVATE USE_GPL)
endif()

if(MINIBWA_ENABLE_MIMALLOC)
	add_library(minibwa_mimalloc OBJECT third_party/mimalloc/static.c)
	target_include_directories(minibwa_mimalloc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/mimalloc)
	target_compile_definitions(minibwa_mimalloc PRIVATE
		NDEBUG
		MI_MALLOC_OVERRIDE
		MI_OSX_INTERPOSE=1
		MI_OSX_ZONE=1
	)
	target_compile_options(minibwa_mimalloc PRIVATE -O3 -Wall -Wextra)
	set_target_properties(minibwa_mimalloc PROPERTIES C_STANDARD 11 C_EXTENSIONS ON)
	target_sources(minibwa_cli PRIVATE $<TARGET_OBJECTS:minibwa_mimalloc>)
else()
	target_compile_definitions(minibwa PRIVATE HAVE_KALLOC)
	target_compile_definitions(minibwa_cli PRIVATE HAVE_KALLOC)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
	check_c_compiler_flag("-msse4.2" MINIBWA_HAS_SSE42)
	check_c_compiler_flag("-mpopcnt" MINIBWA_HAS_POPCNT)
	if(MINIBWA_HAS_SSE42)
		target_compile_options(minibwa PRIVATE -msse4.2)
		target_compile_options(minibwa_cli PRIVATE -msse4.2)
	endif()
	if(MINIBWA_HAS_POPCNT)
		target_compile_options(minibwa PRIVATE -mpopcnt)
		target_compile_options(minibwa_cli PRIVATE -mpopcnt)
	endif()
endif()

if(MINIBWA_ENABLE_OPENMP)
	find_package(OpenMP COMPONENTS C)
	if(OpenMP_C_FOUND)
		target_compile_definitions(minibwa_cli PRIVATE LIBSAIS_OPENMP)
		target_link_libraries(minibwa_cli PRIVATE OpenMP::OpenMP_C)
	endif()
endif()

if(MINIBWA_BUILD_EXAMPLES)
	foreach(example IN ITEMS mbmap-one mbmap-batch)
		string(REPLACE "mbmap-" "ex-" example_source "${example}")
		add_executable(${example} examples/${example_source}.c)
		target_include_directories(${example} PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${CMAKE_CURRENT_SOURCE_DIR}/src
		)
		target_link_libraries(${example} PRIVATE minibwa ZLIB::ZLIB m)
	endforeach()
endif()

install(TARGETS minibwa minibwa_cli
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES include/minibwa.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES docs/man/minibwa.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
