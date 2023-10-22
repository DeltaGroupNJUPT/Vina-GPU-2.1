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

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include "scoring_function.h"
#include "matrix.h"
/*
  结构体precalculate_element传入的数据为unint n和double factor
		fast——double向量，n个向量元素，初始化全为0；
		smooth——结构体向量，有两个成员变量first、second，两个变量都是double类型，n个向量元素，初始化全为0；
		factor——double类型，等于输入的factor。
*/
//理解pair
// pair<int, double> p1;  使用默认构造函数，当一个函数需要返回2个数据的时候，选用pair实现是一个结构体，
//主要的两个成员变量是为first、second 因为是使用struct不是class，所以可以直接使用pair的成员变量
//p1.first = 1;
//p1.second = 2.5;
//cout << p1.first << p1.second << endl;
//输出1  2.5

struct precalculate_element {
	precalculate_element(sz n, fl factor_) : fast(n, 0), smooth(n, pr(0, 0)), factor(factor_) {}

	/*
	input——r2：double
	inter——i：unint
	output——fast[i]：double向量中第i位的数值；
	fast——double向量，向量元素初始化全为0；
	fast.size()计算fast向量元素个数；
	首先先判断r2 * factor是否小于fast向量的元素个数，是的话就进行执行下面的程序，否则终止程序；
	然后对factor * r2强制取整得到i；
	最后判断i 是否小于fast.size()，是的话就进行执行下面的程序即返回fast[i]的值，否则终止程序。
	*/
	fl eval_fast(fl r2) const {
		assert(r2 * factor < fast.size());//它的条件返回错误，则终止程序执行
		sz i = sz(factor * r2);           //r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i < fast.size());          //它的条件返回错误，则终止程序执行
		return fast[i];
	}


	/*
	input——r2：double；
	inter——r2_factored、rem：double， i1、i2：unint, p1、p2：构造函数；
	output——pr(e, dor)：返回两个数值，都为double类型；
	smooth.size()计算fast向量元素个数；
	smooth——结构体向量，有两个成员变量first、second，两个变量都是double类型，n个向量元素，初始化全为0；
	epsilon_fl——计算机能识别的浮点数可表示的最小值，我的计算机能识别的最小浮点数是2.22045e-16

	首先先判断r2_factored=r2 * factor是否小于smooth向量的元素个数，是的话就进行执行下面的程序，否则终止程序；
	然后对r2_factored取整得到i1，而i2=i1+1；
	然后判断i1、i2 是否小于smooth.size()，是的话就进行执行下面的程序即计算rem为factor * r2的小数数据，否则终止程序。
	最后判断rem是否大于-epsilon_fl且小于1 + epsilon_fl，是的话将结构体向量smooth[i1]、smooth[i2]依次引用为p1、p2，所以有p1.first、p1.second、p2.first、p2.second，否则终止程序。
				e   = smooth[sz(factor * r2)].first  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].first  - smooth[sz(factor * r2)].first );
				dor = smooth[sz(factor * r2)].second  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].second   - smooth[sz(factor * r2)].second  );

	*/

	pr eval_deriv(fl r2) const {
		fl r2_factored = factor * r2;         
		assert(r2_factored + 1 < smooth.size());              //它的条件返回错误，则终止程序执行
		sz i1 = sz(r2_factored);
		sz i2 = i1 + 1;                                      // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i1 < smooth.size());
		assert(i2 < smooth.size());
		fl rem = r2_factored - i1;
		assert(rem >= -epsilon_fl);                          //epsilon_fl计算机能识别的浮点数可表示的最小值，我的计算机能识别的最小浮点数是2.22045e-16
		assert(rem < 1 + epsilon_fl);
		const pr& p1 = smooth[i1];	                         //	这里的&不是赋值地址，而是一种引用，也可以说是别名			
		const pr& p2 = smooth[i2];                           //
		fl e   = p1.first  + rem * (p2.first  - p1.first);      
		fl dor = p1.second + rem * (p2.second - p1.second); 
		return pr(e, dor);
	}




/*
input——rs：double向量；
inter——n：unint，dor、delta、r、f1、f2：double
smooth.size()计算smooth向量元素个数；
smooth——结构体向量，有两个成员变量first、second，两个变量都是double类型，向量元素初始化全为0；
fast——double向量，向量元素初始化全为0；
fast.size()计算fast向量元素个数；
首先计算smooth向量元素个数并赋值给n；
然后判断输入rs向量元素个数和fast向量元素个数是否都为n，否的话就终止程序，是的话继续下面程序；
最后通过一个0<=i<n循环，求dor、fast向量元素，其中dor为smooth[i].second，修改dor就等修改smooth[i].second：


		计算dor：首先将向量的第一个元素和最后一个元素都给值为0，否则：

						dor=（smooth[i+1].first - smooth[i-1].first）/（rs[i+1] - rs[i-1]）*rs[i]）；

	    计算fast向量：首先当i == n-1时f2=0；否则f2=smooth[i+1].first，所以：

				     	fast[i]=(smooth[i].first+f2)/2。
*/

	void init_from_smooth_fst(const flv& rs) {
		sz n = smooth.size();
		VINA_CHECK(rs.size() == n);														//assert（）
		VINA_CHECK(fast.size() == n);												    //assert（）
		VINA_FOR(i, n) {																//for(uint VINA_MACROS_TMP=n,i=0; i<VINA_MACROS_TMP;i++)循环
			// calculate dor's
			fl& dor = smooth[i].second;                                                 //这里的&不是赋值地址，而是一种引用，也可以说是别名，修改dor就等修改smooth[i].second	
			if(i == 0 || i == n-1)                                                      //向量的第一个元素和最后一个元素都为0
				dor = 0;
			else {
				fl delta = rs[i+1] - rs[i-1];                                           //求两个间隔一个元素的向量数值的差值，并把这个间隔的元素赋值给 r
				fl r = rs[i];
				dor = (smooth[i+1].first - smooth[i-1].first) / (delta * r);		   
			}
			// calculate fast's from smooth.first's
			fl f1 = smooth[i].first;
			fl f2 = (i+1 >= n) ? 0 : smooth[i+1].first;
			fast[i] = (f2 + f1) / 2;
		}
	}


/*			查找smooth.first中最小值的索引
	input——无；
	inter——i：unint；
	output——tmp：unint；
	smooth.size()计算smooth向量元素个数；
	smooth——结构体向量，有两个成员变量first、second，两个变量都是double类型，向量元素初始化全为0；
	
	先对tmp赋初值0；
	最后通过一个0<=i_inv<smooth.size循环对tmp进行赋值，具体过程如下：

								    i = smooth.size() - i_inv - 1; 

					当i_inv == 0 或是 smooth[i].first < smooth[tmp].first时：

									tmp = i;
*/

	sz min_smooth_fst() const {
		sz tmp = 0;                                                 // returned if smooth.empty()
		VINA_FOR_IN(i_inv, smooth) {                               //for(uint VINA_MACROS_TMP=smooth.size,i_inv=0; i_inv<VINA_MACROS_TMP;i_inv++)
			sz i = smooth.size() - i_inv - 1;                      // i_inv < smooth.size()  => i_inv + 1 <= smooth.size()
			if(i_inv == 0 || smooth[i].first < smooth[tmp].first)
				tmp = i;
		}
		return tmp;
	}


	/*
	input——rs：double向量；left、right：double；
	inter——tmp：与smooth同一元素数量的double向量， min_index：unint， optimal_r、r：double；

	smooth.size()计算smooth向量元素个数；
	smooth——结构体向量，有两个成员变量first、second，两个变量都是double类型，向量元素初始化全为0；

	首先是调用上面min_smooth_fst()查找smooth.first中最小值的索引赋值给min_index，
	并判断输入的double向量rs的元素个数是否大于这个索引值，以及是否smooth的元素个数相等，是的话教训下面程序，否则终止程序；
	然后通过0<=i<smooth.size循环对tmp[i]赋值：
					    //这里rs[i]修改不会修改输入的rs[i]，因为是通过r操作；
						if (rs[i] < rs[min_index] - left ) 
								rs[i] = r+left;
						else if(rs[i] > rs[min_index]+ right)
								rs[i] = rs[i]-right;
						else  
						        rs[i] = rs[min_index];

						if(rs[i] < 0) 
							rs[i] = 0;
						if(rs[i] > rs.back())
						    rs[i] = rs.back();                   //rs.back() 返回rs最末一个元素

						tmp[i] = eval_deriv(rs[i]^2).first;      //eval_deriv()函数返回pr(e, dor)，所以.first应该取值e的值

	最后通过0<=i<smooth.size循环对smooth[i].first]赋值：

						smooth[i].first = tmp[i];

	*/

	void widen_smooth_fst(const flv& rs, fl left, fl right) {
		flv tmp(smooth.size(), 0);								// the new smooth[].first
		sz min_index = min_smooth_fst();						// smooth.first中最小值的索引
		VINA_CHECK(min_index < rs.size());					    //assert // won't hold for n == 0
		VINA_CHECK(rs.size() == smooth.size());					//assert
		fl optimal_r   = rs[min_index];							//
		VINA_FOR_IN(i, smooth) {                                //for(uint VINA_MACROS_TMP=smooth.size,i=0; i<VINA_MACROS_TMP;i++)
			fl r = rs[i];
			if     (r < optimal_r - left ) r += left;
			else if(r > optimal_r + right) r -= right;
			else                           r = optimal_r;

			if(r < 0) r = 0;
			if(r > rs.back()) r = rs.back();                   //back() 返回最末一个元素

			tmp[i] = eval_deriv(sqr(r)).first;                 //eval_deriv()函数返回pr(e, dor)
		}
		VINA_FOR_IN(i, smooth)									//for(uint VINA_MACROS_TMP=smooth.size,i=0; i<VINA_MACROS_TMP;i++)                         
			smooth[i].first = tmp[i];
	}



    /*这个函数主要调用了上面widen_smooth_fst函数和init_from_smooth_fst函数；
	input——rs：double向量，left、right：double；
	*/

	void widen(const flv& rs, fl left, fl right) {
		widen_smooth_fst(rs, left, right);
		init_from_smooth_fst(rs);
	}


	flv fast;                           //double向量
	prv smooth;                         //结构体向量
	fl factor;                          //double
};



/*
      input——sf：scoring_function结构体;   v:计算机能识别的double数可表示的最大值，我的计算机能识别的最大double数是00DB17EE,
	  factor_:double类型初始值32， 
	  data——模板参数为precalculate_element结构体的triangular_matrix结构体；
	  m_cutoff_sqr=weighted_terms结构体中重构的cutoff（）函数的返回值cutoff_的平方；
	  n=unint（m_cutoff_sqr*factor_）+3；factor=factor_(初始值)；
	  data结构体中向量m_data={原子数类型unint*（原子数类型unint-1）/2，输入（n,factor_）为precalculate_element结构体}；
	  m_dim=原子数类型unint；
	  m_atom_typing_used是原子数类型枚举类型；
		最开始要保证factor大于计算机能识别的浮点数可表示的最小值（2.22045e-16），并且n要大于sz(m_cutoff_sqr*factor) + 1和m_cutoff_sqr*factor + 1 
	  

*/
struct precalculate {
	    precalculate(const scoring_function& sf, fl v = max_fl, fl factor_ = 32) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
		m_cutoff_sqr(sqr(sf.cutoff())),     //返回weighted_terms结构体中重构的cutoff（）函数的返回值cutoff_的平方
		n(sz(factor_ * m_cutoff_sqr) + 3),  // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
		factor(factor_),

		data(num_atom_types(sf.atom_typing_used()), precalculate_element(n, factor_)),//num_atom_types（）函数返回weighted_terms结构体中atom_typing_used_所表示原子数类型unint
			//strictly_triangular_matrix(sz n, const T& filler_val) : m_data(n*(n - 1) / 2, filler_val), m_dim(n) {}
			
			m_atom_typing_used(sf.atom_typing_used()) {                                   //

		VINA_CHECK(factor > epsilon_fl);              //factor大于计算机能识别的浮点数可表示的最小值，我的计算机能识别的最小浮点数是2.22045e-16
		VINA_CHECK(sz(m_cutoff_sqr*factor) + 1 < n); // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
		VINA_CHECK(m_cutoff_sqr*factor + 1 < n);

		flv rs = calculate_rs();                                                   //rs[i]=（i/factor）^(1/2);(0-n)
														                          
		VINA_FOR(t1, data.dim())                                                  //for(t1=0;t1<m_dim;t1++)  m_dim=原子数类型unint
			VINA_RANGE(t2, t1, data.dim()) {                                      //for(t2=t1;i<m_dim;t2++) 嵌套循环 m_dim！次
				precalculate_element& p = data(t1, t2);	                         //引用 data（t1, t2）结构体                        
				// init smooth[].first					                          
				VINA_FOR_IN(i, p.smooth)                                          //for(i=0;i<smooth.size;i++)
					p.smooth[i].first = (std::min)(v, sf.eval(t1, t2, rs[i]));   // sf.eval(t1, t2, rs[i]返回乘累加结果，
																				  //并取计算机能取到最大double值和乘累加结果的最小值

				// init the rest
				p.init_from_smooth_fst(rs);                                       //初始化
			}
	}

/*	
	input——type_pair_index：unint，r2：double；
    m_cutoff_sqr=weighted_terms结构体中重构的cutoff（）函数的返回值cutoff_的平方；
	保证输入的r2<= m_cutoff_sqr；
	最后返回fast[sz(factor * r2)},	precalculate_element结构体中fast——double向量，向量元素初始化全为0；
*/
	fl eval_fast(sz type_pair_index, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		return data(type_pair_index).eval_fast(r2);
	}

	/*
		input——type_pair_index：unint，r2：double；
		m_cutoff_sqr=weighted_terms结构体中重构的cutoff（）函数的返回值cutoff_的平方；
		保证输入的r2<= m_cutoff_sqr；
		最后返回(e, dor)；
		e   = smooth[sz(factor * r2)].first  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].first  - smooth[sz(factor * r2)].first );
		dor = smooth[sz(factor * r2)].second  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].second   - smooth[sz(factor * r2)].second  );

	*/

	pr eval_deriv(sz type_pair_index, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		return data(type_pair_index).eval_deriv(r2);
	}


	sz index_permissive(sz t1, sz t2) const { return data.index_permissive(t1, t2); }//返回 (t1 < t2) ? index(t1, t2) : index(t2, t1)
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }             //返回 原子数类型枚举类型
	fl cutoff_sqr() const { return m_cutoff_sqr; }                                   //返回m_cutoff_sqr=weighted_terms结构体中重构的cutoff（）函数的返回值
																				    //cutoff_的平方；

	/*
	
	
	*/
	void widen(fl left, fl right) {
		flv rs = calculate_rs();                           //rs[i]=（i/factor）^(1/2);(0-n)   
		VINA_FOR(t1, data.dim())                          //for(t1=0;t1<m_dim;t1++)；m_dim=原子数类型unint
			VINA_RANGE(t2, t1, data.dim())                //for(t2=t1;i<m_dim;t2++) 嵌套循环 m_dim！次
				data(t1, t2).widen(rs, left, right);     //调用了上面widen_smooth_fst函数和init_from_smooth_fst函数
	}


public:

	/* 
	函数返回double向量tmp，tmp中的值为（i/factor）^(1/2)
	*/
	flv calculate_rs() const {
		flv tmp(n, 0);
		VINA_FOR(i, n)											//for(uint VINA_MACROS_TMP=n,i=0; i<VINA_MACROS_TMP;i++)
			tmp[i] = std::sqrt(i / factor);                     //tmp[i]=(i / factor)^(1/2)
		return tmp;
	}



	fl m_cutoff_sqr;                                          //double
	sz n;													  //unint 
	fl factor;												  //double
	atom_type::t m_atom_typing_used;                          //atom_type::t枚举型

	triangular_matrix<precalculate_element> data;             //模板参数为precalculate_element结构体的triangular_matrix结构体

};

#endif
