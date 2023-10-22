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

#ifndef VINA_TEE_H
#define VINA_TEE_H

#include <iostream>
#include "file.h"


/*


*/
struct tee {
	ofile* of;
	tee() : of(NULL) {}                        //不创建文件

/*
     重载ofile结构体of，用于创建一个名字由name赋值的文件
*/
	void init(const path& name) {
		of = new ofile(name);                 
	}

	/*
	      定义了一个虚函数用于删除of结构体
	*/
	virtual ~tee() { delete of; } 


	/*
	  flush() 是把缓冲区的数据强行输出, 主要用在IO中，即清空缓冲区数据，
	  一般在读写流(stream)的时候，数据是先被读到了内存中，再把数据写到文件中，
	  当你数据读完的时候不代表你的数据已经写完了，因为还有一部分有可能会留在内存这个缓冲区中。
	  这时候如果你调用了close()方法关闭了读写流，那么这部分数据就会丢失，所以应该在关闭读写流之前先flush()。
	  若创建了一个文件将留在缓冲区的数据写到创建的文件中
	*/

	void flush() {
		std::cout << std::flush;                            
		if(of)                      
			(*of) << std::flush;
	}



	/*
	endl() 换行并刷新缓冲区，与flush（）差不多
	*/

	void endl() {                 
		std::cout << std::endl;
		if(of)
			(*of) << std::endl;
	}

	/*
	通过该函数实现of结构体中文件按a独立的控制符号格式输出，符号解释如下链接：
	https://blog.csdn.net/yedawei_1/article/details/105538439
	*/

	void setf(std::ios::fmtflags a) {   //格式控制信息的枚举类型fmtflags ，影响到如何解释输入序列的格式、
		                                //如何生成输出序列的格式，例如整数是16进制还是10进制表示，浮点数是科学计数法还是定点形式等；
										//a用来定义输入输出格式
		std::cout.setf(a);
		if(of)
			of->setf(a);
	}
	/*
	通过该函数实现of结构体中文件按a，b组合控制符号格式输出，符号解释如下链接：
	https://blog.csdn.net/yedawei_1/article/details/105538439
	*/
	void setf(std::ios::fmtflags a, std::ios::fmtflags b) {
		std::cout.setf(a, b);
		if(of)
			of->setf(a, b);
	}
};



/*重构tee结构体
	input ——out：tee结构体，x：任意类型；
	输出x，并且若创建文件就将x写进文件，最后返回out结构体



*/
template<typename T>
tee& operator<<(tee& out, const T& x) {
	std::cout << x;
	if(out.of)
		(*out.of) << x;
	return out;
}

#endif
