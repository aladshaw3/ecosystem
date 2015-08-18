//----------------------------------------
//  Created by Austin Ladshaw on 07/27/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 
 DISCLAIMER: Niether Austin Ladshaw, nor the Georgia Institute of Technology, is the author or owner of any YAML Library or source code.
 Only the files labeld "yaml_tests" were created by Austin Ladshaw for the sole purpose of running and testing the yaml code before
 implementation in the main adsorption software packages developed by Austin Ladshaw at the Georgia Institute of Technology.
 
 The YAML Library (LibYAML) was written by Kirill Simonov and is released under the MIT license. For more information on YAML, go to
 pyyaml.org/wiki/LibYAML. The MIT License is provided below...
 
 */

/*
 
 The MIT License (MIT)
 
 Copyright (c) 2015 Austin Ladshaw
 Portions copyright 2006 Kirill Simonov 
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 
 */

#include "yaml_wrapper.h"

#ifndef YAML_TESTS_HPP_
#define YAML_TESTS_HPP_

typedef enum
{
	DOCUMENT_HEADER,
	SUB_HEADER,
	KEY
	
} yaml_key_type;

typedef struct
{
	int key_type = -1;
	int key_type_old = -1;
	int value_type = -1;
	int value_type_old = -1;
	
} yaml_read_data;

int YAML_TEST01();

int YAML_TEST02();

int YAML_TEST03();

#endif
