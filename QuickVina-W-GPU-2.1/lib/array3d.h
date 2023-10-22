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

#ifndef VINA_ARRAY3D_H
#define VINA_ARRAY3D_H

#include <exception> // std::bad_alloc
#include "common.h"

/*
intput:		uint i,j
function:	进行i*j=tmp运算,结合判断返回值tmp。返回类型为uint
			a.输入中有0则输出0;
			b.结果溢出抛出异常;
			c.正常输出
*/
inline sz checked_multiply(sz i, sz j) {
	if(i == 0 || j == 0) return 0;				//判断运算中是否有0以减少运算量?
	const sz tmp = i * j;						//计算i*j
	if(tmp < i || tmp < j || tmp / i != j)		//溢出判断
		throw std::bad_alloc(); // can't alloc if the size makes sz wrap around
	return tmp;
}

/* 
intput:		uint i,j,k
function:   计算不溢出且输入非零的i*j*k 
解释：嵌套 checked_multiply 并返回 uint型结果
*/
inline sz checked_multiply(sz i, sz j, sz k) {
	return checked_multiply(checked_multiply(i, j), k);
}



template<typename T>							//C++模板,形参(不需要明确定义变量类型)
/*
侵入式(intrusive)序列化
参考:https://blog.csdn.net/zj510/article/details/8105408
*/
class array3d {		
public:
	sz m_i, m_j, m_k;							//uint
	std::vector<T> m_data;						//T类型名为m_data空对象
	friend class boost::serialization::access;
	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & m_i;
		ar & m_j;
		ar & m_k;
		ar & m_data;
	}
public:
	array3d() : m_i(0), m_j(0), m_k(0) {}
	array3d(sz i, sz j, sz k) : m_i(i), m_j(j), m_k(k), m_data(checked_multiply(i, j, k)) {}
	sz dim0() const { return m_i; }
	sz dim1() const { return m_j; }
	sz dim2() const { return m_k; }
	sz dim(sz i) const {
		switch(i) {
			case 0: return m_i;
			case 1: return m_j;
			case 2: return m_k;
			default: assert(false); return 0; // to get rid of the warning
		}
	}
	void resize(sz i, sz j, sz k) { // data is essentially garbled
		m_i = i;
		m_j = j;
		m_k = k;
		m_data.resize(checked_multiply(i, j, k));
	}
	T&       operator()(sz i, sz j, sz k)       { return m_data[i + m_i*(j + m_j*k)]; }
	const T& operator()(sz i, sz j, sz k) const { return m_data[i + m_i*(j + m_j*k)]; }
};

#endif
