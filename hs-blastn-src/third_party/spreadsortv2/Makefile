spreadsort: example/sample.cpp 
		g++ -Wno-deprecated -O3 example/sample.cpp -o spreadsort -Wall -I .; chmod 777 spreadsort

alreadysorted: example/alreadysorted.cpp 
		g++ -Wno-deprecated -O3 example/alreadysorted.cpp -o alreadysorted -Wall -I .; chmod 777 alreadysorted

mostlysorted: example/mostlysorted.cpp 
		g++ -Wno-deprecated -O3 example/mostlysorted.cpp -o mostlysorted -Wall -I .; chmod 777 mostlysorted

rightshift: example/rightshiftsample.cpp 
		g++ -Wno-deprecated -O3 example/rightshiftsample.cpp -o rightshift -Wall -I .; chmod 777 rightshift

reverseintsort: example/reverseintsample.cpp 
		g++ -Wno-deprecated -O3 example/reverseintsample.cpp -o reverseintsort -Wall -I .; chmod 777 reverseintsort

int64: example/int64.cpp 
		g++ -Wno-deprecated -O3 example/int64.cpp -o int64 -Wall -I .; chmod 777 int64

floatsort: example/floatsample.cpp 
		g++ -Wno-deprecated -O3 example/floatsample.cpp -o floatsort -Wall -I .; chmod 777 floatsort

shiftfloatsort: example/shiftfloatsample.cpp 
		g++ -Wno-deprecated -O3 example/shiftfloatsample.cpp -o shiftfloatsort -Wall -I .; chmod 777 shiftfloatsort

floatfunctorsort: example/floatfunctorsample.cpp 
		g++ -Wno-deprecated -O3 example/floatfunctorsample.cpp -o floatfunctorsort -Wall -I .; chmod 777 floatfunctorsort

double: example/double.cpp 
		g++ -Wno-deprecated -O3 example/double.cpp -o double -Wall -I .; chmod 777 double

stringsort: example/stringsample.cpp 
	g++ -Wno-deprecated -O3 example/stringsample.cpp -o stringsort -Wall -I .; chmod 777 stringsort

wstringsort: example/wstringsort.cpp 
		g++ -Wno-deprecated -O3 example/wstringsortsample.cpp -o wstringsort -Wall -I .; chmod 777 wstringsort

reversestringsort: example/reversestringsample.cpp  
	g++ -Wno-deprecated -O3 example/reversestringsample.cpp -o reversestringsort -Wall -I .; chmod 777 reversestringsort

charstringsort: example/charstringsample.cpp  
	g++ -Wno-deprecated -O3 example/charstringsample.cpp -o charstringsort -Wall -I .; chmod 777 charstringsort

stringfunctorsort: example/stringfunctorsample.cpp  
	g++ -Wno-deprecated -O3 example/stringfunctorsample.cpp -o stringfunctorsort -Wall -I .; chmod 777 stringfunctorsort

reversestringfunctorsort: example/reversestringfunctorsample.cpp  
	g++ -Wno-deprecated -O3 example/reversestringfunctorsample.cpp -o reversestringfunctorsort -Wall -I .; chmod 777 reversestringfunctorsort

keyplusdata: example/keyplusdatasample.cpp 
	g++ -Wno-deprecated -O3 example/keyplusdatasample.cpp -o keyplusdata -Wall -I .; chmod 777 keyplusdata
	
randomgen: example/randomGen.cpp
	g++ -Wno-deprecated -O3 example/randomGen.cpp -o randomgen; chmod 777 randomgen

all: randomgen spreadsort rightshift floatsort stringsort stringfunctorsort keyplusdata floatfunctorsort charstringsort reversestringsort reversestringfunctorsort shiftfloatsort reverseintsort alreadysorted mostlysorted int64 double

clean:
		rm -f *.o spreadsort spreadsort.exe firstdim firstdim.exe randomgen randomgen.exe keyplusdata keyplusdata.exe floatsort floatsort.exe stringfunctorsort stringfunctorsort.exe rightshift rightshift.exe charstringsort charstringsort.exe stringsort stringsort.exe floatfunctorsort floatfunctorsort.exe reversestringsort reversestringsort.exe reversestringfunctorsort reversestringfunctorsort.exe shiftfloatsort shiftfloatsort.exe reverseintsort reverseintsort.exe alreadysorted alreadysorted.exe mostlysorted mostlysorted.exe int64 int64.exe double double.exe