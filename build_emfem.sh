#/bin/bash

# EMFEM 编译脚本
# 适用于通过Spack安装依赖并手动下载nanoflann的情况

# 1. 进入源码目录（请将 /path/to/emfem/sources 替换为您的实际路径）
SOURCE_DIR="/home/dt/桌面/AMS1/ams_2"
cd "$SOURCE_DIR" || { echo "错误：无法进入源码目录 $SOURCE_DIR"; exit 1; }

# 2. 创建并进入构建目录
mkdir -p build
cd build || { echo "错误：无法进入build目录"; exit 1; }

# 3. 加载Spack环境（基于第一张图，您的spack在/home/dt/spack）
export PATH=$HOME/spack/bin:$PATH
source $HOME/spack/share/spack/setup-env.sh

# 4. 加载所需的Spack包
spack load cmake
spack load mpi
spack load metis
spack load tetgen
spack load eigen
spack load petsc

# 5. 手动下载nanoflann（如果尚未下载）
NANOFLANN_FILE="$HOME/nanoflann.hpp"
if [ ! -f "$NANOFLANN_FILE" ]; then
    echo "正在下载nanoflann.hpp..."
    wget https://raw.githubusercontent.com/jlblancoc/nanoflann/master/include/nanoflann.hpp -O "$NANOFLANN_FILE"
    if [ $? -eq 0 ]; then
        echo "nanoflann.hpp 下载成功"
    else
        echo "警告：nanoflann下载失败，请手动下载或检查网络"
        echo "请手动执行: wget https://raw.githubusercontent.com/jlblancoc/nanoflann/master/include/nanoflann.hpp -O ~/nanoflann.hpp"
        exit 1
    fi
fi

# 6. 运行CMake配置（基于第二张图格式）
echo "正在配置CMake..."
cmake -DCMAKE_BUILD_TYPE=Release \
    -DPETSC_DIR=$(spack location -i petsc) \
    -DTETGEN_DIR=$(spack location -i tetgen) \
    -DNANOFLANN_DIR=$(dirname "$NANOFLANN_FILE") \
    -DEIGEN_DIR=$(spack location -i eigen) \
    -DMETIS_DIR=$(spack location -i metis) ../src_ams 

# 7. 编译项目
echo "正在编译..."
make -j$(nproc)

echo "编译完成！"
