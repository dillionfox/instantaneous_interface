xdrlibdir="/home/dillion/pkg/surf-master/depends/lib"
xdrincldir="/home/dillion/pkg/surf-master/depends/xdrfile-1.1.4/include"
/opt/intel/bin/icpc -DCPLUSPLUS -I $xdrincldir -L $xdrlibdir -lxdrfile -o chandler_willard.out chandler_willard.cpp -Wall
