cmake -S . -G Ninja -B build \
    -DSKBUILD_PROJECT_NAME="diptest" \
    -DSKBUILD_PROJECT_VERSION="0.8.0" \
    -DDIPTEST_MBUILD=ON \
    -DDIPTEST_CPP_STANDARD="11" \
    -DDIPTEST_ENABLE_DEVMODE=ON \
    -DDIPTEST_ENABLE_DEBUG=OFF \
    -DDIPTEST_ENABLE_OPENMP=ON \
    -DDIPTEST_ENABLE_EXT_TESTS=OFF \
    -DDIPTEST_ENABLE_ARCH_FLAGS=ON \
    -DDIPTEST_ENABLE_64BIT_INDEX=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DOpenMP_ROOT=$(brew --prefix)/opt/libomp \
    -Dpybind11_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())") \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

cmake --build build --target install --config Release --parallel 4
