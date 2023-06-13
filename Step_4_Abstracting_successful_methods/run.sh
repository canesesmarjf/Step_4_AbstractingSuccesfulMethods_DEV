if [ $(uname) = "Linux" ]; then
  export LD_LIBRARY_PATH=/home/jfcm/Repos/binaryTree/arma_libs/lib:$LD_LIBRARY_PATH
fi

FORGE_PATH=/home/jfcm/arm/forge/22.0.1/bin
export PATH=${FORGE_PATH}:$PATH

if [ "$1" = "2" ]; then
  ./bin/main_2.exe
  echo "Completed main_2.exe"
elif [ "$1" = "3" ]; then
  ./bin/main_3.exe
  echo "Completed main_3.exe"
else
  echo "Invalid argument. Usage: $0 [2|3]"
  exit 1
fi
