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

#ifndef VINA_PARALLEL_PROGRESS_H
#define VINA_PARALLEL_PROGRESS_H

#include <boost/progress.hpp>
#include <boost/thread/mutex.hpp>

#include "incrementable.h"


/*parallel_progress类，继承incrementable
* 私有成员mutex类self（互斥锁），progress_display类指针p
* 构造函数parallel_progress初始化p为空指针
* 成员函数init：
*       input:unsigned long数n
*       output:无
*       作用：申明进度条，可以在控制台上显示程序的执行进度，参数即为进度条一行的总个数n
*		      该进度条保存在堆上，p指向进度条的地址，在局部函数调用结束后仍然能够使用
* 成员函数operator，重载++运算符：
*        input:无
*        output:无
*        作用：1.区域锁将互斥锁封装到对象self_lk，避免多个线程同时访问p 【可能解释有误】
*              2.计算++parallel_progress时，如果p不为空指针，并将p指针所指的进度条自增1
* 析构函数parallel_progress（可以被子类继承和覆盖）：
*        input:无
*        output:无
*        作用：释放new分配的单个对象p指针指向的内存
*/

struct parallel_progress : public incrementable {
	parallel_progress() : p(NULL) {}
	void init(unsigned long n) { p = new boost::progress_display(n); }
	void operator++() {
		if(p) {
			boost::mutex::scoped_lock self_lk(self);
			++(*p);
		}
	}
	virtual ~parallel_progress() { delete p; }
private:
	boost::mutex self;
	boost::progress_display* p;
};

#endif

