/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_STATISTICS_H
#define VINA_STATISTICS_H

#include <algorithm> // sort
#include <cmath> // sqrt
#include "common.h" // for flv

// simple but inefficient implementations


/*   求均值
	input—— v:double向量, ,只读；
	inter——acc:double;
	output——mean:double;
	v.size()计算v向量元素个数；
	v.empty()用于判断v向量是否为空，为空返回1否则为0；
	        acc=v[i]+acc;i=0-v.size
			if(v向量为空)
			mean=0;
			else 
			mean=acc/v.size;
*/

inline fl mean(const flv& v) {
	fl acc = 0;
	VINA_FOR_IN(i, v)						  //for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += v[i];						  //对向量中的值进行累加
	return v.empty() ? 0 : (acc/v.size());    //若v为空向量返回0，否则返回累加的值除向量元素个数
}




/*   求偏差
input—— v:double向量, ,只读；
inter——acc:double、m:double;
output——deviation:double;
v.size()计算v向量元素个数；
v.empty()用于判断v向量是否为空，为空返回1否则为0；
m = mean(v)用于求v向量的均值
             acc=(v[i]-mean(m))^2+acc;i=0-v.size
             if(v向量为空)
             deviation=0;
             else
             deviation=(acc/v.size)^(1/2);
*/

inline fl deviation(const flv& v) {
	fl m = mean(v);
	fl acc = 0;
	VINA_FOR_IN(i, v)								//for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += sqr(v[i] - m);
	return v.empty() ? 0 : std::sqrt(acc/v.size());
}


/*   求均方差根
input—— a、b:double向量 ,只读；
inter——acc:double、m:double;
output——rmsd:double;
a.size()计算a向量元素个数；
a.empty()用于判断a向量是否为空，为空返回1否则为0；
计算前要先判断a、b向量元素个数是否相等，相等才能进行接下来的运算。否则退出；
              acc=(a[i]-b[i])^2+acc;i=0-a.size
              if(v向量为空)
              rmsd=0;
              else
              rmsd=(acc/a.size)^(1/2);
*/

inline fl rmsd(const flv& a, const flv& b) {
	VINA_CHECK(a.size() == b.size());                 //断言 assert(a.size() == b.size())，若a、b向量个数相同继续执行，否则终止程序
	fl acc = 0;
	VINA_FOR_IN(i, a)                                //for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += sqr(a[i] - b[i]);
	return a.empty() ? 0 : std::sqrt(acc/a.size());
}



/*  
input—— a、b:double向量 ,只读；
inter——acc:double、m:double;
output——average_difference:double;
a.size()计算a向量元素个数；
a.empty()用于判断a向量是否为空，为空返回1否则为0；
计算前要先判断a、b向量元素个数是否相等，相等才能进行接下来的运算。否则退出；
           acc=(b[i]-a[i])+acc;i=0-a.size
           if(v向量为空)
           average_difference=0;
           else
           average_difference=(acc/a.size);
*/

inline fl average_difference(const flv& b, const flv& a) { // b - a
	VINA_CHECK(a.size() == b.size());                      //断言 assert(a.size() == b.size())，若a、b向量个数相同继续执行，否则终止程序
	fl acc = 0;
	VINA_FOR_IN(i, a)									   //for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += b[i] - a[i];
	return a.empty() ? 0 : (acc/a.size());
}





/*   
input—— x、y:double向量,只读；
inter——sum_x,sum_y,sum_x_sq,sum_y_sq,sum_prod、sd_x、sd_y、cov、tmp:double;
output——pearson:double;
x.size()计算a向量元素个数；
x.empty()用于判断a向量是否为空，为空返回1否则为0；
epsilon_fl为计算机能识别的浮点数可表示的最小值，2.22045e-16；
计算前要先判断a、b向量元素个数是否相等，相等才能进行接下来的运算。否则退出；
若向量为空直接返回0；
                     sum_x=sum_x+x[i];i=0-n(a.size)
					 sum_y=sum_y+y[i];
					 sum_x_sq=sum_x_sq+x[i]^2;
					 sum_y_sq=sum_y_sq+y[i]^2;
					 sum_prod=sum_prod+y[i]*x[i];
					 sd_x=(sum_x_sq/n - (sum_x/n)^2)^(1/2);//标准差
					 sd_y=(sum_y_sq/n - (sum_y/n)^2)^(1/2);//标准差
					 cov= sum_prod/n - (sum_x/n) * (sum_y/n);//协方差
					 tmp=sd_x * sd_y;
					 if（abs(tmp) < epsilon_fl）
					 pearson= 0;
					 else 
					 pearson= cov / tmp;

*/

inline fl pearson(const flv& x, const flv& y) {
	sz n = x.size();
	VINA_CHECK(n == y.size());                          //断言 assert(a.size() == b.size())，若a、b向量个数相同继续执行，否则终止程序
	if(n == 0) return 0; 
	fl sum_x = 0;
	fl sum_y = 0;
	fl sum_x_sq = 0;
	fl sum_y_sq = 0;
	fl sum_prod = 0;
	
	VINA_FOR(i, n) {                                  //for(uint VINA_MACROS_TMP=n,i=0; i<VINA_MACROS_TMP;i++)
		sum_x += x[i];
		sum_y += y[i];
		sum_x_sq += sqr(x[i]);
		sum_y_sq += sqr(y[i]);
		sum_prod += x[i] * y[i];
	}
	fl sd_x = std::sqrt(sum_x_sq/n - sqr(sum_x/n));  // FIXME the argument is supposed to be >= 0, but ...
	fl sd_y = std::sqrt(sum_y_sq/n - sqr(sum_y/n));
	fl cov = sum_prod/n - (sum_x/n) * (sum_y/n);
	fl tmp = sd_x * sd_y;
	if(std::abs(tmp) < epsilon_fl) return 0;        // epsilon_fl为计算机能识别的浮点数可表示的最小值，2.22045e-16
	return cov / tmp;
}





/*
x:double;
i:unsigned int;
将传到spearman_aux中的值传给x、i；
*/
struct spearman_aux {
	fl x;
	sz i;
	spearman_aux(fl x, sz i) : x(x), i(i) {}
};




/*这个函数为下面的get_rankings排序函数提供排序方式——升序
将函数中的形参改为指针a、b，可以减少内存空间，
加入const之后，在这个operator函数中一旦对a.x、b.x、a.i、b.i有修改的操作就会报错，
可以防止我们的误操作
	if(a.x<b.x)
	return 1；
	else
	return 0;
*/
inline bool operator<(const spearman_aux& a, const spearman_aux& b) {     //
	return a.x < b.x;
}


/*实现对输入x向量中的值进行升序排序，并记录其坐标，将在坐标存进tmp向量中

input—— x:double向量,只读；
output——tmp:double向量,为x.size，初始化全为0的向量;
inter——to_sort：结构体向量daxiao ；
首先to_sort结构体x变量中在最后面依次添加输入x[i]（i=0-x.size），即将结构体中的x赋值输入向量x;
to_sort结构体i变量中在最后面依次添加输入i（i=0-x.size），即对结构体中的i赋值;
然后结合上一个operator符号重载函数，对x进行升序；
最后在按照x升序顺序取i的值，使得tmp[i]=rank :rank=0-x.size;

*/

inline flv get_rankings(const flv& x) {
	std::vector<spearman_aux> to_sort;                 //to_sort为结构体向量
	VINA_FOR_IN(i, x)								   //for(uint VINA_MACROS_TMP=x.size,i=0; i<VINA_MACROS_TMP;i++)
	to_sort.push_back(spearman_aux(x[i], i));      // push_back函数:将一个新的元素加到vector的最后面，位置为当前最后一个元素的下一个元素
	std::sort(to_sort.begin(), to_sort.end());         //（1）第一个参数first：是要排序的数组的起始地址。
													   //（2）第二个参数last：是结束的地址（最后一个数据的后一个数据的地址）
													   //（3）第三个参数comp是排序的方法：可以是从升序也可是降序。
													   //如果第三个参数不写，则默认的排序方法是从小到大排序。

	flv tmp(x.size(), 0);							   //定义一个大小为x.size，初始化全为0的向量
	VINA_FOR_IN(rank, to_sort)						   //for(uint VINA_MACROS_TMP=to_sort.size,rank=0; rank<VINA_MACROS_TMP;rank++)
		tmp[to_sort[rank].i] = rank;
	return tmp;
}

/*
input—— x、y:double向量,只读；
output——pearson：double；
x、y进行升序后返回各自的坐标向量（get_rankings作用）x1，y1(这里为了区分所以自己定义的)；
最后对向量x1，y1进行数学运算（pearson作用），返回结果。
*/

inline fl spearman(const flv& x, const flv& y) {
	return pearson(get_rankings(x), get_rankings(y));
}

#endif
