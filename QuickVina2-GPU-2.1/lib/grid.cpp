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

#include "grid.h"
/*input ： 三个结构体组成的结构体组
function：给结构体中的五个数组赋值
output：void    改变数组数值

*/
void grid::init(const grid_dims& gd) {//输入为三个结构体组成的结构体组
	m_data.resize(gd[0].n+1, gd[1].n+1, gd[2].n+1);   //三个结构体中的项n+1然后乘   ？？？？？？？？？
	m_init = vec(gd[0].begin, gd[1].begin, gd[2].begin);//给类grid中的数组赋值     每个结构体的begin
	m_range = vec(gd[0].span(), gd[1].span(), gd[2].span());//给类grid中的数组赋值  每个结构体的end - begin
	assert(m_range[0] > 0);//确保数组数字>0
	assert(m_range[1] > 0);
	assert(m_range[2] > 0);
	m_dim_fl_minus_1 = vec(m_data.dim0() - 1.0, 
	                       m_data.dim1() - 1.0,
			               m_data.dim2() - 1.0);//数组三个数都减一   给类grid中的数组赋值
	VINA_FOR(i, 3) {
		m_factor[i] = m_dim_fl_minus_1[i] / m_range[i];   //m_dim_fl_minus_1/m_range    给类grid中的数组赋值
		m_factor_inv[i] = 1 / m_factor[i];//1/m_factor     给类grid中的数组赋值
	}
}// 
/*input ：长度3double类型数组location  只读     double：  slope、v   指针类型deriv/NULL
output：  double
inter：长度3double类型数组s   miss  长度3int类型数组region  长度3uint类型数组a
function：输出double类型的数并且改变deriv对应的值  见文件    //分为指针有效和无效类型  

*/
fl grid::evaluate_aux(const vec& location, fl slope, fl v, vec* deriv) const { // sets *deriv if not NULL
	vec s  = elementwise_product(location - m_init, m_factor); //a0*b0    a1*b1  a2*b2   构造数组1   根据类中的数组和输入给数组s赋值

	vec miss(0, 0, 0);    //构造数组2
	boost::array<int, 3> region;  //构造数组3
	boost::array<sz, 3> a;  //构造数组4

	VINA_FOR(i, 3) {
		if(s[i] < 0) {    //根据数组1   给数组1234赋值
			miss[i] = -s[i];
			region[i] = -1;
			a[i] = 0; 
			s[i] = 0;
		}       
		else if(s[i] >= m_dim_fl_minus_1[i]) {//根据数组1   给数组1234赋值
			miss[i] = s[i] - m_dim_fl_minus_1[i];
			region[i] = 1;
			assert(m_data.dim(i) >= 2);
			a[i] = m_data.dim(i) -  2; 
			s[i] = 1;
		}
		else {
			region[i] = 0; // now that region is boost::array, it's not initialized
			a[i] = sz(s[i]); //根据数组1   给数组1234赋值
			s[i] -= a[i];
		}
		assert(s[i] >= 0);
		assert(s[i] <= 1);
		assert(a[i] >= 0);
		assert(a[i]+1 < m_data.dim(i));//判断s[i]在0-1之间  a[i]>=0且a[0]<m_i a[1]<m_j a[2]<m_k
	}
	const fl penalty = slope * (miss * m_factor_inv); // 输出的一部分

	const sz x0 = a[0];
	const sz y0 = a[1];
	const sz z0 = a[2];

	const sz x1 = x0+1;
	const sz y1 = y0+1;
	const sz z1 = z0+1;


	const fl f000 = m_data(x0, y0, z0);
	const fl f100 = m_data(x1, y0, z0);
	const fl f010 = m_data(x0, y1, z0);
	const fl f110 = m_data(x1, y1, z0);
	const fl f001 = m_data(x0, y0, z1);
	const fl f101 = m_data(x1, y0, z1);
	const fl f011 = m_data(x0, y1, z1);
	const fl f111 = m_data(x1, y1, z1);

	const fl x = s[0];
	const fl y = s[1];
	const fl z = s[2];

	const fl mx = 1-x;
	const fl my = 1-y;
	const fl mz = 1-z;

	fl f = 
		f000 *  mx * my * mz  +
		f100 *   x * my * mz  +
		f010 *  mx *  y * mz  + 
		f110 *   x *  y * mz  +
		f001 *  mx * my *  z  +
		f101 *   x * my *  z  +
		f011 *  mx *  y *  z  +
		f111 *   x *  y *  z  ;

	if(deriv) { // valid pointer  有指针类型参数
		const fl x_g = 
			f000 * (-1)* my * mz  +
			f100 *   1 * my * mz  +
			f010 * (-1)*  y * mz  + 
			f110 *   1 *  y * mz  +
			f001 * (-1)* my *  z  +
			f101 *   1 * my *  z  +
			f011 * (-1)*  y *  z  +
			f111 *   1 *  y *  z  ;


		const fl y_g = 
			f000 *  mx *(-1)* mz  +
			f100 *   x *(-1)* mz  +
			f010 *  mx *  1 * mz  + 
			f110 *   x *  1 * mz  +
			f001 *  mx *(-1)*  z  +
			f101 *   x *(-1)*  z  +
			f011 *  mx *  1 *  z  +
			f111 *   x *  1 *  z  ;


		const fl z_g =  
			f000 *  mx * my *(-1) +
			f100 *   x * my *(-1) +
			f010 *  mx *  y *(-1) + 
			f110 *   x *  y *(-1) +
			f001 *  mx * my *  1  +
			f101 *   x * my *  1  +
			f011 *  mx *  y *  1  +
			f111 *   x *  y *  1  ;

		vec gradient(x_g, y_g, z_g);
		curl(f, gradient, v);//改变f和gradient
		vec gradient_everywhere;

		VINA_FOR(i, 3) {//改变deriv对应的值
			gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
			(*deriv)[i] = m_factor[i] * gradient_everywhere[i] + slope * region[i];
		}

		return f + penalty;//返回输出
	}
	else {  //无输入指针类型参数
		curl(f, v);
		return f + penalty;
	}
} 
