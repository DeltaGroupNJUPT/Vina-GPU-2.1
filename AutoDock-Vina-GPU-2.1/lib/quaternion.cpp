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

#include "quaternion.h"
/*	判断是否四元数标准化
                        input——q：四元数主要用于旋转，
                        返回值为bool类型的函数，
                        eq(boost::math::abs(q), 1)   ——abs(boost::math::abs(q) - 1)<0.001，返回1否则0，abs(q) 为系数的平方相加开根号，这里使用的是系统函数计算四元数绝对值;
                        eq(quaternion_norm_sqr(q), 1)——abs(quaternion_norm_sqr(q) - 1)<0.001，返回1否则0，quaternion_norm_sqr(q)计算四元数系数的平方相加结果
                        两者进行相与返回最终结果
*/
bool quaternion_is_normalized(const qt& q) { // not in the interface, used in assertions
	return eq(quaternion_norm_sqr(q), 1) && eq(boost::math::abs(q), 1);
}


/*	判断四元数是否近似相等
	input——a、b：四元数；

	假设const qt qt_identity(1, 2, 3, 4);
	R_component_1()//四元数第一个元素1;
	R_component_2()//四元数第一个元素2;
	R_component_3()//四元数第一个元素3;
	R_component_4()//四元数第一个元素4;

	eq(a.R_component_1(), b.R_component_1())这个函数如果a.R_component_1()-b.R_component_1()的差的绝对值小于0.001则返回true，否则false


*/

bool eq(const qt& a, const qt& b) { // elementwise approximate equality - may return false for equivalent rotations
	return eq(a.R_component_1(), b.R_component_1()) && \
		eq(a.R_component_2(), b.R_component_2()) && \
		eq(a.R_component_3(), b.R_component_3()) && \
		eq(a.R_component_4(), b.R_component_4());
}



/*	角度四元数
	input——axis：vec结构体；angle：double；
	首先判断axis结构体中(data[0]^2+data[1]^2+data[2]^2)^(1/2)与1的差的绝对值是否小于0.001，是的话继续，否则终止（是否标准化）
	normalize_angle(angle)用来角度归一化；
	返回四元数（cos（angle/2），sin（angle/2）*data[0]，sin（angle/2）*data[1]，sin（angle/2）*data[2]）
*/
qt angle_to_quaternion(const vec& axis, fl angle) { // axis is assumed to be a unit vector
													//assert(eq(tvmet::norm2(axis), 1));
	assert(eq(axis.norm(), 1));                    //判断axis结构体中(data[0]^2+data[1]^2+data[2]^2)^(1/2)与1的差的绝对值是否小于0.001，是的话继续，否则终止
	normalize_angle(angle); // this is probably only necessary if angles can be very big //角度归一化，只在角度很大的时候才可能有用
	fl c = std::cos(angle / 2);                    //cos（x/2）
	fl s = std::sin(angle / 2);						//sin（x/2）
	return qt(c, s*axis[0], s*axis[1], s*axis[2]);  //返回四元数，axis[0]=data[0],axis[1]=data[1],axis[2]=data[2]
}



/*   pre_angle_to_quaternion
	input——rotation：vec结构体；
	inter——angle：double；axis——vec结构体；
	output——四元数；
	angle=(data[0]^2+data[1]^2+data[2]^2)^(1/2)；data数组是vec结构体成员变量
		if（angle>2.22045e-16）
		重构的四元数乘法运算axis = ((1/angle)*data[0](1/angle)*data[1],(1/angle)*data[2])；
		然后通过angle_to_quaternion(axis, angle)函数返回角度四元数
		else
		返回实数1
*/
qt angle_to_quaternion(const vec& rotation) {
	//fl angle = tvmet::norm2(rotation); 
	fl angle = rotation.norm();
	if (angle > epsilon_fl) {
		//vec axis; 
		//axis = rotation / angle;	
		vec axis = (1 / angle) * rotation;   //重构的四元数乘法运算((1/angle)*data[0](1/angle)*data[1],(1/angle)*data[2])
		return angle_to_quaternion(axis, angle);
	}
	return qt_identity;
}



/*
	input——q：四元数；
	inter——c、angle、s：double；axis——vec结构体；

	axis结构体中的data数组储存q四元数三个虚数值
	首先判断q是否标准化，然后将q实数赋值给c；
		if（c在(-1, 1)内）
			angle=2*acos(c); //反余弦acos
			if（angle>pi）
				angle=angle-2*pi；
			if（abs（sin(angle / 2)）<2.22045e-16）
				返回 vec(0, 0, 0)；
			data[0]=data[0]*(angle/sin(angle / 2));
			data[1]=data[1]*(angle/sin(angle / 2));
			data[2]=data[2]*(angle/sin(angle / 2));
			返回axis	
	     else
			 返回 vec(0, 0, 0)；
*/

vec quaternion_to_angle(const qt& q) {
	assert(quaternion_is_normalized(q));		//判断是否四元数标准化
	const fl c = q.R_component_1();			    //q实数部分
	if (c > -1 && c < 1) {						// c may in theory be outside [-1, 1] even with approximately normalized q, due to rounding errors
		fl angle = 2 * std::acos(c);		    // acos is in [0, pi]  反余弦
		if (angle > pi)
			angle -= 2 * pi;					// now angle is in [-pi, pi]
		vec axis(q.R_component_2(), q.R_component_3(), q.R_component_4());  //q的虚数部分
		fl s = std::sin(angle / 2);				 // perhaps not very efficient to calculate sin of acos
		if (std::abs(s) < epsilon_fl)
			return zero_vec;
		axis *= (angle / s);                     //重构*
		return axis;
	}
	else // when c = -1 or 1, angle/2 = 0 or pi, therefore angle = 0
		return zero_vec;
}



/*
	input——q：四元数；
	inter——a、b、c、d、aa、ab、ac、ad、bb、bc、bd、cc、cd、dd：double；
	output——tmp：mat结构体；
	首先判断q是否四元数标准化，然后a等于q的实数，b等于第一个虚数值，c等于第二个虚数值，d等于第三个虚数值；
	aa = a*a;ab = a*b;ac = a*c;ad = a*d;bb = b*b
	bc = b*c;bd = b*d;cc = c*c;cd = c*d;dd = d*d
	接着判断aa + bb + cc + dd是否近似等于1，是的话继续，否则终止；
	mat结构体中：
	data[0]=(aa + bb - cc - dd)；data[3]=2 * (-ad + bc)     ；data[6]=2 * (ac + bd)      ;
	data[1]=( 2 * (ad + bc)    ；data[4]=(aa - bb + cc - dd)；data[7]= 2 * (-ab + cd)    ;
	data[2]=2 * (-ac + bd)     ；data[5]=2 * (ab + cd)      ；data[8]=(aa - bb - cc + dd);
*/

mat quaternion_to_r3(const qt& q) {
	assert(quaternion_is_normalized(q));   //判断q是否四元数标准化

	const fl a = q.R_component_1();
	const fl b = q.R_component_2();
	const fl c = q.R_component_3();
	const fl d = q.R_component_4();

	const fl aa = a*a;
	const fl ab = a*b;
	const fl ac = a*c;
	const fl ad = a*d;
	const fl bb = b*b;
	const fl bc = b*c;
	const fl bd = b*d;
	const fl cc = c*c;
	const fl cd = c*d;
	const fl dd = d*d;

	assert(eq(aa + bb + cc + dd, 1));

	mat tmp;

	// from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
	tmp(0, 0) = (aa + bb - cc - dd);
	tmp(0, 1) = 2 * (-ad + bc);
	tmp(0, 2) = 2 * (ac + bd);

	tmp(1, 0) = 2 * (ad + bc);
	tmp(1, 1) = (aa - bb + cc - dd);
	tmp(1, 2) = 2 * (-ab + cd);

	tmp(2, 0) = 2 * (-ac + bd);
	tmp(2, 1) = 2 * (ab + cd);
	tmp(2, 2) = (aa - bb - cc + dd);

	return tmp;
}


/*随机方向
	input——generator：随机数生成器
	inter——q：四元数；nrm：double；
	random_normal(0, 1, generator)，在double类型的输入0、1之间生成正态分布的随机数，并返回函数
	q的四个元素为随机数；
	nrm为q的绝对值，假设q（a,b,c,d）即qd的绝对值为(a^2+b^2+c^2+d^2)^(1/2)
	if（nrm>2.22045e-16)
	q(a/nrm,b/nrm,c/nrm,d/nrm);若近似标准化返回q
	else
	继续生成随机方向，这一步应该不会发生
*/
qt random_orientation(rng& generator) {
	qt q(random_normal(0, 1, generator),
		random_normal(0, 1, generator),
		random_normal(0, 1, generator),
		random_normal(0, 1, generator));
	fl nrm = boost::math::abs(q);      //算出四元数的绝对值即系数平方相加开根号
	if (nrm > epsilon_fl) {           //epsilon_fl=2.22045e-16
		q /= nrm;
		assert(quaternion_is_normalized(q));
		return q;
	}
	else
		return random_orientation(generator); // this call should almost never happen这个不应该发生
}


/*
	input——q：四元数；rotation——vec结构体
	首先判断q是否标准化，
	是的话q等于angle_to_quaternion(rotation)返回的四元数和q的乘积
	四元数乘法：q1（w1,x1,y1,z1）;q2（w2,x2,y2,z2）
	q1 * q2 =
	(w1*w2 - x1*x2 - y1*y2 - z1*z2) +
	(w1*x2 + x1*w2 + y1*z2 - z1*y2) i +
	(w1*y2 - x1*z2 + y1*w2 + z1*x2) j +
	(w1*z2 + x1*y2 - y1*x2 + z1*w2) k

	最后将q近似标准化
*/
void quaternion_increment(qt& q, const vec& rotation) {
	assert(quaternion_is_normalized(q));              //是否标准化
	q = angle_to_quaternion(rotation) * q;
	quaternion_normalize_approx(q); // normalization added in 1.1.2
									//quaternion_normalize(q); // normalization added in 1.1.2
}


/*
	input——a、b：四元数;
	inter——tmp：四元数；
	依次判断a、b是否标准化，然后将b复制给tmp；
	tmp=tmp/a;
	四元数除法：q1/q2——q1（w1,x1,y1,z1）;q2（w2,x2,y2,z2）
	 deno = w2*w2+x2*x2+y2*y2+z2*z2；
 ((w1*w2+x1*x2+y1*y2+z1*z2)/deno,(-w1*x2+x1*w2-y1*z2+z1*y2)/deno,(-w1*y2+x1*z2+y1*w2-z1*x2)/deno,(-w1*z2-x1*y2+y1*x2+z1*w2)/deno)
   最后返回quaternion_to_angle(tmp)；
*/
vec quaternion_difference(const qt& b, const qt& a) { // rotation that needs to be applied to convert a to b
	quaternion_is_normalized(a);
	quaternion_is_normalized(b);
	qt tmp = b;
	tmp /= a; // b = tmp * a    =>   b * inv(a) = tmp 
	return quaternion_to_angle(tmp); // already assert normalization
}


/*
	打印出四元数，类似于qt(a,b,c,d)

*/
void print(const qt& q, std::ostream& out) { // print as an angle
	print(quaternion_to_angle(q), out);
}
