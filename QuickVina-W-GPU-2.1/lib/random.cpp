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

#include <ctime> // for time (for seeding)

#include "random.h"
#include "my_pid.h"



/*input  a、b   ：double
	out：random_fl   ：double
function:	在double类型的输入a、b之间生成均匀分布的随机数，并返回函数  double类型
			a.保证a在B的左边
			b.调用bosst的函数  一个是用来生成随机数  另一个是制定区间   最后构造随机数生成器
			c.验证随机数并返回double类型的随机数
*/


fl random_fl(fl a, fl b, rng& generator) { // rng用于生成 随机数
	assert(a < b); // 保证区间
	typedef boost::uniform_real<fl> distr;// 指定均匀分布的区间，包含两个端点，取double类型
	boost::variate_generator<rng&, distr> r(generator, distr(a, b));  // 构造符合要求的随机数生成器
	fl tmp = r();
	assert(tmp >= a); //验证生成随机数在区间内
	assert(tmp <= b);
	return tmp;
}



/*input  mean、sigma   ：double
	out：random_normal   ：double
function:	在double类型的输入mean、sigma之间生成正态分布的随机数，并返回函数  double类型
			a.保证sigma不小于0
			b.调用bosst的函数  一个是用来生成随机数  另一个是制定符合要求区间   最后构造随机数生成器
			c.返回double类型的随机数
*/

fl random_normal(fl mean, fl sigma, rng& generator) { 
	assert(sigma >= 0); // 正态分布保证sigma大于等于0
	typedef boost::normal_distribution<fl> distr; //定义了一个默认返回 double 型浮点值的正态分布
	boost::variate_generator<rng&, distr> r(generator, distr(mean, sigma));// 构造符合要求的随机数生成器
	return r();
}


/*input  a、b   ：int
	out：random_int   ：int
function:	在int类型的输入a、b之间生成均匀分布的随机整数，并返回函数  int类型
			a.保证a在B的左边
			b.调用bosst的函数  一个是用来生成随机数  另一个是制定区间   最后构造随机数生成器
			c.验证随机数并返回int类型的随机数
*/

//  产生一个在ab之间的随机数  整数类型
int random_int(int a, int b, rng& generator) { // rng用于生成 随机数
	assert(a <= b); // 保证区间
	typedef boost::uniform_int<int> distr;  // 指定均匀分布的区间，包含两个端点，取整数值
	boost::variate_generator<rng&, distr> r(generator, distr(a, b));   // 构造符合要求的随机数生成器
	int tmp = r();
	assert(tmp >= a);    //验证生成随机数在区间内
	assert(tmp <= b);
	return tmp;
}
/*input  a、b   ：unsigned int
	out：random_int   ：unsigned int
function:	在unsigned int类型的输入a、b之间生成均匀分布的随机整数，并返回函数 unsigned int类型
			a.保证a在B的左边  并且a、b不小于0
			b.调用bosst的函数  一个是用来生成随机数  另一个是制定区间   最后构造随机数生成器
			c.验证随机数并返回int类型的随机数
*/
sz random_sz(sz a, sz b, rng& generator) { // 产生一个在ab之间的随机数  类型
	assert(a <= b);
	assert(int(a) >= 0);
	assert(int(b) >= 0);   //保证a在B的左边  并且a、b不小于0
	int i = random_int(int(a), int(b), generator);    //生成随机数
	assert(i >= 0);
	assert(i >= int(a));
	assert(i <= int(b));    // 保证随机数在区间内
	return static_cast<sz>(i);
}
/*input  true   ： int
	out：tmp数组   ：double
function:	当信号ture到来的时候，返回一个数组，并且数组的三个数在以0为圆心，半径为1的球体内  double
			a.生成三个随机数  double
			b.组成数组
			c.验证在球内
*/
vec random_inside_sphere(rng& generator) {
	while(true) { // on average, this will have to be run about twice
		fl r1 = random_fl(-1, 1, generator);    //产生-1——1的随机浮点数
		fl r2 = random_fl(-1, 1, generator);
		fl r3 = random_fl(-1, 1, generator);

		vec tmp(r1, r2, r3);//把这三个数组成数组
		if(sqr(tmp) < 1)//r1的平方+r2的平方+r3的平方再开根号<1  确保在球内
			return tmp;
	}
}
/*input   corner1、corner2    ： double类型  数组
	out：tmp数组   ：double
function:	输入两个数组，产生一个随机数组，要求生成在corner1[i]-corner2[i]之间   double

*/
vec random_in_box(const vec& corner1, const vec& corner2, rng& generator) { // expects corner1[i] < corner2[i]
	vec tmp;  //定义数组
	VINA_FOR_IN(i, tmp)     //for循环  i从0到数组的长度
		tmp[i] = random_fl(corner1[i], corner2[i], generator);//生成随机数组
	return tmp;
}
/*
	out：auto_seed   ：unsigned int
function:	返回时间和当前进程的标示符(PID)的乘积   int类型？？？？
*/
int auto_seed() { 
	return my_pid() * int(std::time(NULL));
}
