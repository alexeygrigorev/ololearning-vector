include_directories(include)

file(GLOB SOURCES 
    "src/*.cpp"
    "src/la/math_utils.h"
    "src/la/math_utils.cpp"
    "src/la/DenseVector.cpp"
    "src/la/DenseMatrix.cpp"
    "src/la/svd.cpp"
    )

# aux_source_directory("src" SOURCE)

add_library(libololearning SHARED STATIC ${SOURCES})
target_include_directories(libololearning PUBLIC include)
