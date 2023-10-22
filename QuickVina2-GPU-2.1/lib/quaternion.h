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

#ifndef VINA_QUATERNION_H
#define VINA_QUATERNION_H

#include <boost/math/quaternion.hpp>
#include <boost/serialization/split_free.hpp>

#include "common.h"
#include "random.h"

typedef boost::math::quaternion<fl> qt;              //四元数
bool quaternion_is_normalized(const qt & q);         //返回bool类型的quaternion_is_normalized函数，输入函数为四元数

// non-intrusive free function split serialization


/*
		save函数：input——qt:四元数
		将四元数的一个实数三个虚数存储在文件中或者通过网络发出去。

		load函数：将文件中的四元数的一个实数三个虚数写进qt四元数。
		序列化调用save函数，反序列化调用load函数。
		假设const qt qt_identity(1, 2, 3, 4);
		R_component_1()//四元数第一个元素1;
		R_component_2()//四元数第一个元素2;
		R_component_3()//四元数第一个元素3;
		R_component_4()//四元数第一个元素4;


*/
namespace boost {
	namespace serialization {
		template<class Archive>
		void save(Archive& ar, const qt& q, const unsigned version) {//假设const qt qt_identity(1, 2, 3, 4);
			fl q1 = q.R_component_1();						//q1:四元数第一个元素1;
			fl q2 = q.R_component_2();						//q2:四元数第一个元素2;
			fl q3 = q.R_component_3();						//q3:四元数第一个元素3;
			fl q4 = q.R_component_4();						//q4:四元数第一个元素4;

			ar & q1;
			ar & q2;
			ar & q3;
			ar & q4;
		}


		template<typename Archive>
		void load(Archive& ar, qt& q, const unsigned version) {
			fl a, b, c, d;
			ar & a;
			ar & b;
			ar & c;
			ar & d;
			q = qt(a, b, c, d);
		}
	}
}


BOOST_SERIALIZATION_SPLIT_FREE(qt)//BOOST_SERIALIZATION_SPLIT_FREE()来告诉Boost序列化库使用save和load代替serialize函数




/*
	各种声明
*/
  
bool eq(const qt& a, const qt& b); // elementwise approximate equality - may return false for equivalent rotations
const qt qt_identity(1, 0, 0, 0);   //定义了一个实部为1，虚部都为0的四元数
qt angle_to_quaternion(const vec& axis, fl angle); // axis is assumed to be a unit vector
qt angle_to_quaternion(const vec& rotation); // rotation == angle * axis
vec quaternion_to_angle(const qt& q);
mat quaternion_to_r3(const qt& q);

/*    四元数绝对值的平方

		input——qt:四元数
		若qt=a+bi+cj+dk;
		函数返回：a^2+b^2+c^2+d^2

*/
inline fl quaternion_norm_sqr(const qt& q) { // equivalent to sqr(boost::math::abs(const qt&))
	return sqr(q.R_component_1()) + sqr(q.R_component_2()) + sqr(q.R_component_3()) + sqr(q.R_component_4());
}



/*	四元数标准化

	input——q：四元数
	这个函数主要通过quaternion_norm_sqr（）函数求得的四元数系数的平方相加后在将其开根号求四元数绝对值
	若求得的值大于计算机能识别的浮点数可表示的最小值2.22045e-16，则：
		q=q*（1/a）
	最后判断bool(abs(boost::math::abs(q) - 1)<0.001)与上bool(abs(quaternion_norm_sqr(q) - 1)<0.001)，为真继续程序否则终止程序，
*/

inline void quaternion_normalize(qt& q) {
	const fl s = quaternion_norm_sqr(q);
	assert(eq(s, sqr(boost::math::abs(q))));  //如果自己写的系数平方和与系统的函数计算的四元数绝对值平方的差的绝对值小于0.001则返回true，否则false   
    const fl a = std::sqrt(s);                //		求四元数绝对值
	assert(a > epsilon_fl);                   //epsilon_fl计算机能识别的浮点数可表示的最小值，我的计算机能识别的最小浮点数是2.22045e-16
	q *= 1/a;                                 //
	assert(quaternion_is_normalized(q));
}


/*四元数近似标准化
		input——q：四元数，tolerance：double，1e-6；
		若系数的平方和减1的绝对值小于1e-6，则无操作，保持初始值，否则判断系数的平方和的开跟号是否大于计算机能识别的浮点数可表示的最小值2.22045e-16，是的话：
			q=q*（1/a）
	最后判断bool(abs(boost::math::abs(q) - 1)<0.001)与上bool(abs(quaternion_norm_sqr(q) - 1)<0.001)，为真继续程序否则终止程序，


*/
inline void quaternion_normalize_approx(qt& q, const fl tolerance = 1e-6) {
	const fl s = quaternion_norm_sqr(q);      //系数的平方和
	assert(eq(s, sqr(boost::math::abs(q))));  //如果自己写的系数平方和与系统的函数计算的四元数绝对值平方的差的绝对值小于0.001则返回true，否则false   
    if(std::abs(s - 1) < tolerance)            //系数的平方和减一的绝对值小于1e-6
        ; // most likely scenario
    else {
        const fl a = std::sqrt(s);             //求四元数绝对值
        assert(a > epsilon_fl);
        q *= 1/a;
        assert(quaternion_is_normalized(q));
    }
}

/*
各种声明
*/

qt random_orientation(rng& generator);
void quaternion_increment(qt& q, const vec& rotation);
vec quaternion_difference(const qt& b, const qt& a); // rotation that needs to be applied to convert a to b
void print(const qt& q, std::ostream& out = std::cout); // print as an angle

#endif
