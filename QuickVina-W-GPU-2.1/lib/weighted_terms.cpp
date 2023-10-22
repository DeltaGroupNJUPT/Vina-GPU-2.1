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

#include "weighted_terms.h"


/*
	input——t：terms结构体、weights：double向量；double类型cutoff_初始化为0，输入传入weighted_terms结构体中；
	inter——enabled_usable_terms：unint向量 ，atom_typing_used_：atom_type::t枚举类型；
		首先判断模板参数为distance_additive的term_set结构体中enabled中大于0的值的个数、模板参数为additive的term_set结构体中enabled中大于0的值的个数
	和模板参数为intermolecular的term_set结构体中enabled中大于0的值的个数 是否都等于0；是的话继续程序，否则退出；
		if（模板参数为usable的term_set结构体中enabled[i]大于0）(i=0-模板参数为usable的term_set结构体中fun向量元素个数）
			if（enabled_usable_terms向量为空）{
				atom_typing_used_等于usable_terms.size为模板参数为usable中fun[i]的atom_typing_used枚举型；
		    else
			    判断atom_typing_used_是否等于usable_terms.size为模板参数为usable中fun[i]的atom_typing_used枚举型；是的话继续否则终止程序；
			在enabled_usable_terms向量后面添加i序号；
			cutoff_赋值为模板参数为usable中fun[i]的cutoff
			}
*/
weighted_terms::weighted_terms(const terms* t, const flv& weights) : t(t), weights(weights), cutoff_(0) { // does not own t
	VINA_CHECK(t->distance_additive_terms.num_enabled() == 0); //该函数返回模板参数为distance_additive的term_set结构体中enabled中大于0的值的个数
	VINA_CHECK(t->         additive_terms.num_enabled() == 0); //该函数返回模板参数为additive的term_set结构体中enabled中大于0的值的个数
	VINA_CHECK(t->   intermolecular_terms.num_enabled() == 0); //该函数返回模板参数为intermolecular的term_set结构体中enabled中大于0的值的个数      
	VINA_FOR_IN(i, t->usable_terms)                            //for(i=0;i<usable_terms.size;i++)，usable_terms.size为模板参数为usable的term_set结构体中fun向量元素个数
		if(t->usable_terms.enabled[i]) {
			if(enabled_usable_terms.empty())					//enabled_usable_terms为unint向量      
				atom_typing_used_ = t->usable_terms[i].atom_typing_used;//atom_typing_used_：atom_type::t枚举类型，usable_terms[i]=fun[i]
			else
				VINA_CHECK(atom_typing_used_ == t->usable_terms[i].atom_typing_used);

			enabled_usable_terms.push_back(i);
			cutoff_ = (std::max)(cutoff_, t->usable_terms[i].cutoff);
		}
}

/*
	input——t1、t2:unint,   r:double;
	inter——acc：double，enabled_usable_terms：unint向量，weights：double向量；
	for(i=0;i<enabled_usable_terms.size;i++)
   对weight[i]和 t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r) 的乘积进行累加
   这里有一个问题t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r)是返回0？？？？？？？？？？？？？？？？
*/

fl weighted_terms::eval(sz t1, sz t2, fl r) const { // intentionally not checking for cutoff
	fl acc = 0;
	VINA_FOR_IN(i, enabled_usable_terms)                  //for(i=0;i<enabled_usable_terms.size;i++)
		acc += weights[i] * t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r);     // t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r)返回0？？？？
	return acc;
}




/*
	input——m:model结构体,e:double;
	inter——enabled_usable_terms：unint向量，weights：double向量，it：迭代器，in：参数为model结构体的conf_independent_inputs结构体；
	output——tmp：dopuble
	it迭代器指向weights向量中第enabled_usable_terms向量元素个数的元素；
	tmp等于e + read_iterator(i) * in.******   // read_iterator(i)为迭代器指向下一个元素 ，******应该是取决于调用这个函数的前一个everthing.cpp中的函数决定的
	最后判断it是否指向weights向量的最后一位，是的话返回tmp
*/

fl weighted_terms::conf_independent(const model& m, fl e) const {
	flv::const_iterator it = weights.begin() + enabled_usable_terms.size();     //指向weights向量中第enabled_usable_terms向量元素个数的元素
	conf_independent_inputs in(m); // FIXME quite inefficient, but I think speed is irrelevant here, right?
	fl tmp = t->eval_conf_independent(in, e, it);
	assert(it == weights.end());        //判断it是否指向weights向量的最后一位
	return tmp;
}
