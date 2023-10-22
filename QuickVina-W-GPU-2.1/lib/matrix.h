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

#ifndef VINA_MATRIX_H
#define VINA_MATRIX_H

#include <vector>
#include "triangular_matrix_index.h"

// these 4 lines are used 3 times verbatim - defining a temp macro to ease the pain  运算符重载
#define VINA_MATRIX_DEFINE_OPERATORS \
	const T& operator()(sz i) const { return m_data[i]; } \
	      T& operator()(sz i)       { return m_data[i]; } \
	const T& operator()(sz i, sz j) const { return m_data[index(i, j)]; } \
	      T& operator()(sz i, sz j)       { return m_data[index(i, j)]; } 
/*
定义结构体matrix 
声明向量m_data  类型为模板参数T
定义uint型函数index   输入  uint  i  j  输出  uint  i + m_i*j
定义构造函数并赋值


定义函数resize   input  uint  m  n    模板参数T类型  filler_val 只读  
function  ：改变m_data  m_i m_j的值
判断
循环执行  ij从00 到 m_im_j   把参数为（i  j）的对象赋值给 tmp[i+m*j]
对m_data  m_i m_j  赋值


定义函数append  input  类matrix类型模板参数为T的参数x   只读      模板参数T类型的参数filler_val  输出viod
function  对参数为(i+m, j+n)的对象赋值

定义uint类型的函数dim_1dim_2  function 返回 m_i  m_j
*/
template<typename T>
class matrix {
	std::vector<T> m_data;//向量 m_data
	sz m_i, m_j;
public:
	sz index(sz i, sz j) const {//input   uint  i j   output  uint  i + m_i*j
		assert(j < m_j);
		assert(i < m_i); 
		return i + m_i*j; // column-major
	}
	matrix() : m_i(0), m_j(0) {}
	matrix(sz i, sz j, const T& filler_val) : m_data(i*j, filler_val), m_i(i), m_j(j) {}
	void resize(sz m, sz n, const T& filler_val) { // new sizes should be the same or greater than the old
		if(m == dim_1() && n == dim_2()) return; // 判断m n   
		VINA_CHECK(m >= dim_1());
		VINA_CHECK(n >= dim_2());
		std::vector<T> tmp(m*n, filler_val);//M*N个元素  每个都是filler_val
		VINA_FOR(i, m_i)
			VINA_FOR(j, m_j)
				tmp[i+m*j] = (*this)(i, j);//  把参数为（i  j）的对象赋值给 tmp[i+m*j]
		m_data = tmp;			//赋值m_data 
		m_i = m;
		m_j = n;
	}
	void append(const matrix<T>& x, const T& filler_val) {
		sz m = dim_1();
		sz n = dim_2();
		resize(m + x.dim_1(), n + x.dim_2(), filler_val);//对m_data  m_i m_j  赋值
		VINA_FOR(i, x.dim_1())
			VINA_FOR(j, x.dim_2())(i+m, j+n)
				(*this)(i+m, j+n) = x(i, j);//对参数为(i+m, j+n)的对象赋值
	}
	VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above
	sz dim_1() const { return m_i; }
	sz dim_2() const { return m_j; }
};
/*定义类triangular_matrix    模板参数为T
* 
函数index  uint类型  input  uint i  j  return  i + j*(j+1)/2;只读
函数index_permissive  uint类型   input  uint  i  j   只读  返回i + j*(j+1)/2;or//j + i*(i+1)/2;
函数dim   uint类型  返回  uint类型   m_dim
*/
template<typename T>
class triangular_matrix {
public:
	std::vector<T> m_data;//向量  T类型
	sz m_dim;
public:
	sz index(sz i, sz j) const { return triangular_matrix_index(m_dim, i, j); }//i + j*(j+1)/2; 
	sz index_permissive(sz i, sz j) const { return (i < j) ? index(i, j) : index(j, i); }//i + j*(j+1)/2;or//j + i*(i+1)/2;
	triangular_matrix() : m_dim(0) {}
	triangular_matrix(sz n, const T& filler_val) : m_data(n*(n+1)/2, filler_val), m_dim(n) {} 
	VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above     宏
	sz dim() const { return m_dim; }
};
/*定义类strictly_triangular_matrix  模板参数 T
*
函数index    uint类型   input  i  j  uint      只读   判断ij关系然后 return   i + j*(j-1)/2

函数index_permissive   uint类型   input  i  j  uint   只读    return   i + j*(j-1)/2 orreturn   j + i*(i-1)/2

函数resize  inout   uint  n    只读T类型filler_val  function  ：给m_dim赋值    将m_data的容量改成*(n-1)/2，扩容的部分设置为filler_val

函数append 	
input m  strictly_triangular_matrix类     filler_val  T类型
out   void
function ：  给类中的参数  m_dim赋值   并且    m_data的容量改成*(n-1)/2，扩容的部分设置为filler_val
	               然后 通过输入更改当前对象的参数  

函数append  	
input   rectangular  matrix类型  只读          triangular     strictly_triangular_matrix类型  只读
out    void
function：      1、判断rectangular的dim_1dim_2，   dim_2==0  直接结束   dim_1=0   把triangular赋值当前对象  结束
				2、执行函数		append(triangular, filler_val)   *给类中的参数  m_dim赋值   并且    当前m_data的容量改成*(n-1)/2，扩容的部分设置为filler_val
	                                             并通过输入triangular更改当前对象的m_data     其中		filler_val=m_data[0]
	            3、然后  输入rectangular的m_data更改当前对象的m_data


*/
template<typename T>
class strictly_triangular_matrix {
	std::vector<T> m_data;  //向量   类型为T
	sz m_dim;
public:
	sz index(sz i, sz j) const {
		assert(j < m_dim);
		assert(i < j); 
		assert(j >= 1);
		return i + j*(j-1)/2;
	}
	sz index_permissive(sz i, sz j) const { return (i < j) ? index(i, j) : index(j, i); }
	strictly_triangular_matrix() : m_dim(0) {}
	strictly_triangular_matrix(sz n, const T& filler_val) : m_data(n*(n-1)/2, filler_val), m_dim(n) {}
	void resize(sz n, const T& filler_val) {
		if(n == m_dim) return; 
		VINA_CHECK(n > m_dim);
		m_dim = n;
		m_data.resize(n*(n-1)/2, filler_val); // preserves original data//将m_data的容量改成*(n-1)/2，扩容的部分设置为filler_val
	}


	void append(const strictly_triangular_matrix<T>& m, const T& filler_val) { 
		sz n = dim();//m_dim
		resize(n + m.dim(), filler_val);//给m_dim赋值    将m_data的容量改成*(n-1)/2，扩容的部分设置为filler_val
		VINA_FOR(i, m.dim())
			VINA_RANGE(j, i+1, m.dim())
				(*this)(i+n, j+n) = m(i, j);//通过输入更改当前对象的m_data  
	}//        m_data[ i+n + m_i*（i+n）]  = m_data[ i + m_i*j]



	void append(const matrix<T>& rectangular, const strictly_triangular_matrix<T>& triangular) {
		VINA_CHECK(dim() == rectangular.dim_1());
		VINA_CHECK(rectangular.dim_2() == triangular.dim());
		if(rectangular.dim_2() == 0) return;//m_i  ==0?
		if(rectangular.dim_1() == 0) {//m_j==0?
			(*this) = triangular;
			return;
		}
		const T& filler_val = rectangular(0, 0); // needed by 'append below'

		sz n = dim();//n=m_dim
		append(triangular, filler_val);
		/*给类中的参数  m_dim赋值   并且    m_data的容量改成*(n-1)/2，扩容的部分设置为filler_val
	               然后 通过输入更改当前对象的m_data*/
		VINA_FOR(i, rectangular.dim_1())   //i=0 ;i<m_i;i++
			VINA_FOR(j, rectangular.dim_2())//j=0 ;j<m_j;j++
				(*this)(i, n + j) = rectangular(i, j);//通过输入rectangular的m_data更改当前对象的m_data
		//   m_data[ i + m_i*（J+n）]  = m_data[ i + m_i*j]
	}
	VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above运算符重载
	sz dim() const { return m_dim; }
};

#undef VINA_MATRIX_DEFINE_OPERATORS

#endif
