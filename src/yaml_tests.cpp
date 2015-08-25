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

#include "yaml_tests.h"
#include <iostream>

//This example shows token based yaml parsing 
int YAML_TEST01()
{
	int success = 0;
	
	//This is the yaml file we want to read in c
	FILE *fh = fopen("public.yaml", "r");
	
	//This is the primary yaml object used by the parser
	yaml_parser_t parser;
	
	//We will first use token based parsing for this example
	//For token based parsing, we need to declare the yaml_token_t structure
	yaml_token_t token;
	
	// Initialize parser and check the file
	//success = yaml_parser_initialize(&parser);
	//if (success == 1) {success = 0;}
	//else {std::cout << "Failed to initialize parser!\n";success = 1;}
	
	//Alternative method
	if(!yaml_parser_initialize(&parser)) //NOTE: this function will return a 1 for success and a 0 for failure
		fputs("Failed to initialize parser!\n", stderr);
	if(fh == NULL)
		fputs("Failed to open file!\n", stderr);
	
	// Set input file - This function will actually open the file. Make sure file exists before opening
	yaml_parser_set_input_file(&parser, fh);
	
	bool key = false;
	bool value = false;
	std::string current_key;
	
	// Loop through the document, read every token, print out the reads
	do {
		
		//Scan through all token that were read one by one in order
		yaml_parser_scan(&parser, &token);
		
		//Switch case that is used to change code behavior based on token type
		switch(token.type)
		{
				/* Stream start/end */
			case YAML_STREAM_START_TOKEN: puts("STREAM START"); key=false; value=false; break;
			case YAML_STREAM_END_TOKEN:   puts("STREAM END");  key=false; value=false; break;
				
				/* Token types (read before actual token) */
			case YAML_KEY_TOKEN:
			{
				printf("(Key token)   ");
				//std::cout << token.data.scalar.value << std::endl; //This is null at this point
				key = true;
				value = false;
				break;
			}
			case YAML_VALUE_TOKEN:
			{
				printf("(Value token) ");
				key = false;
				value = true;
				break;
			}
				
				/* Block delimeters */
			case YAML_BLOCK_SEQUENCE_START_TOKEN: puts("<b>Start Block (Sequence)</b>");key=false; value=false; break;
			case YAML_BLOCK_ENTRY_TOKEN:          puts("<b>Start Block (Entry)</b>");  key=false; value=false;  break;
			case YAML_BLOCK_END_TOKEN:            puts("<b>End block</b>");     key=false; value=false;         break;
				
				/* Data */
			case YAML_BLOCK_MAPPING_START_TOKEN:  puts("[Block mapping]");    key=false; value=false;        break;
			case YAML_SCALAR_TOKEN:
			{
				
				if (key == true) {std::cout << "\n READ IN KEY \n"; char *test2 = (char *)token.data.scalar.value; current_key = test2;}
				if (value == true) {std::cout << "\n READ IN VALUE OF KEY " << current_key << " \n";}
				printf("scalar %s \n", token.data.scalar.value);
				
				double a = 1;
				//a = *reinterpret_cast<double*>(event.data.scalar.value);
				//a = strtod(event.data.scalar.value, NULL);
				//a = *(double*)event.data.scalar.value;
				//a = event.data.scalar.value;
				char *test = (char *) token.data.scalar.value;
				a = atof(test); //use atoi() for ints and will have to check for true and false
				
				int i=0;
				//Function converts array of char to all lower case
				while (test[i])
				{
					test[i] = tolower(test[i]);
					i++;
				}
				std::string word = test;
				
				std::cout << "\ta = " << a << std::endl;
				std::cout << "\tword = " << word << std::endl;
				
				break;
			}
				
				/* Others */
			default:
			{
				printf("Got token of type %d\n", token.type);key=false; value=false;
			}
		}
		
		//This line deletes each token after it has been used
		if(token.type != YAML_STREAM_END_TOKEN)
			yaml_token_delete(&token);
		
	} while(token.type != YAML_STREAM_END_TOKEN);
	yaml_token_delete(&token);
	
	// Cleanup - make sure the file gets closed and the data structure is deleted
	yaml_parser_delete(&parser);
	fclose(fh);
	
	return success;
}

//This example shows event based yaml parsing
int YAML_TEST02()
{
	int success = 0;
	
	//This is the yaml file we want to read in c
	FILE *fh = fopen("public.yaml", "r");
	
	//This is the primary yaml object used by the parser
	yaml_parser_t parser;
	
	//Now we will try event based parsing, which is more "natural"
	//For event based parsing, we need to declare the yaml_event_t structure
	yaml_event_t event;
	
	//Alternative method
	if(!yaml_parser_initialize(&parser)) //NOTE: this function will return a 1 for success and a 0 for failure
		fputs("Failed to initialize parser!\n", stderr);
	if(fh == NULL)
		fputs("Failed to open file!\n", stderr);
	
	// Set input file - This function will actually open the file. Make sure file exists before opening
	yaml_parser_set_input_file(&parser, fh);
	
	// Loop through the document, read every event, print out the reads
	do {
		
		//First, check the event for errors
		if (!yaml_parser_parse(&parser, &event)) {
			printf("Parser error %d\n", parser.error);
			exit(EXIT_FAILURE);
		}
		
		//Switch case for all event types (this is the complete list of events)
		switch(event.type)
		{
				//A default case for when no event occurs
			case YAML_NO_EVENT: puts("No event!"); break;
				
				/* Stream start/end */
			case YAML_STREAM_START_EVENT: puts("STREAM START"); break;
			case YAML_STREAM_END_EVENT:   puts("STREAM END");   break;
				
				/* Block delimeters */
			case YAML_DOCUMENT_START_EVENT: puts("<b>Start Document</b>"); break;
			case YAML_DOCUMENT_END_EVENT:   puts("<b>End Document</b>");   break;
			case YAML_SEQUENCE_START_EVENT: puts("<b>Start Sequence</b>"); break;
			case YAML_SEQUENCE_END_EVENT:   puts("<b>End Sequence</b>");   break;
			case YAML_MAPPING_START_EVENT:  puts("<b>Start Mapping</b>");  break;
			case YAML_MAPPING_END_EVENT:    puts("<b>End Mapping</b>");    break;
				
				/* Data */
			case YAML_ALIAS_EVENT:
				{
					printf("Got alias (anchor %s)\n", event.data.alias.anchor);
					//std::cout << "Value = " << event.data.scalar.value << std::endl; //Can't do this
					//std::cout << "Value = " << event.data.scalar.length << std::endl;
					//event.data.scalar.anchor = event.data.alias.anchor;
					//std::cout << event.data.scalar.value << std::endl;
					break;
				}
			case YAML_SCALAR_EVENT:
			{
				printf("Got scalar (value %s)\n", event.data.scalar.value);
				printf("Got scalar (TAG %s)\n", event.data.scalar.tag);
				double a = 1;
				//a = *reinterpret_cast<double*>(event.data.scalar.value);
				//a = strtod(event.data.scalar.value, NULL);
				//a = *(double*)event.data.scalar.value;
				//a = event.data.scalar.value;
				char *test = (char *) event.data.scalar.value;
				a = atof(test); //use atoi() for ints and will have to check for true and false
				
				int i=0;
				//Function converts array of char to all lower case
				while (test[i])
				{
					test[i] = tolower(test[i]);
					i++;
				}
				std::string word = test;
				
				std::cout << "\ta = " << a << std::endl;
				std::cout << "\tword = " << word << std::endl;
				break;
			}
				
		}
		
		//Delete the current event as you go (optional)
		if(event.type != YAML_STREAM_END_EVENT)
			yaml_event_delete(&event);
		
	} while(event.type != YAML_STREAM_END_EVENT); //Continue loop until last event reached
	
	//Clear out the events after finished using
	yaml_event_delete(&event);
	
	// Cleanup - make sure the file gets closed and the data structure is deleted
	yaml_parser_delete(&parser);
	fclose(fh);
	
	return success;
}

int YAML_TEST03()
{
	int success = 0;
	
	//Declarations
	FILE *fh;
	yaml_parser_t parser;
	yaml_token_t token;
	yaml_token_t token_old;
	yaml_read_data read;
	std::string current_key;
	std::string current_value;
	std::string current_block;
	std::string current_alias;
	std::string current_anchor;
	double d_value;
	bool new_block;
	bool block_key;
	bool sub_block;
	bool doc_start;
	
	//Initializations
	fh = fopen("test_input.yml", "r");
	if(!yaml_parser_initialize(&parser)) //NOTE: this function will return a 1 for success and a 0 for failure
		fputs("Failed to initialize parser!\n", stderr);
	if(fh == NULL)
		fputs("Failed to open file!\n", stderr);
	yaml_parser_set_input_file(&parser, fh);
	new_block = false;
	block_key = false;
	sub_block = false;
	doc_start = false;
	
	// Loop through the document, read every token, print out the reads
	do
	{
		
		//Scan through all token that were read one by one in order
		yaml_parser_scan(&parser, &token);
		
		//Switch case that is used to change code behavior based on token type
		switch(token.type)
		{
			//Cases based on the type of token read in from input
			case YAML_STREAM_START_TOKEN:
			{
				//printf("STREAM START\n");
				break;//Breaks are required
			}
			case YAML_STREAM_END_TOKEN:
			{
				//printf("STREAM END\n");
				break;
			}
			case YAML_KEY_TOKEN:
			{
				printf("\tYAML_KEY_TOKEN\n");
				
				if (token_old.type == YAML_BLOCK_MAPPING_START_TOKEN)
				{
					//About to read a key for the start of a main block
					block_key = true;
				}
				else
				{
					block_key = false;
				}
				
				break;
			}
			case YAML_VALUE_TOKEN:
			{
				printf("\tYAML_VALUE_TOKEN\n");
				break;
			}
			case YAML_BLOCK_SEQUENCE_START_TOKEN:
			{
				//printf("YAML_BLOCK_SEQUENCE_START_TOKEN\n");
				break;
			}
			case YAML_BLOCK_ENTRY_TOKEN:
			{
				printf("YAML_BLOCK_ENTRY_TOKEN\n");
				if (new_block == true)
					sub_block = true;
				else
					sub_block = false;
				new_block = true;
				break;
			}
			case YAML_BLOCK_END_TOKEN:
			{
				printf("YAML_BLOCK_END_TOKEN\n\n");
				if (sub_block == true)
					sub_block = false;
				else
					new_block = false;
				break;
			}
			case YAML_BLOCK_MAPPING_START_TOKEN:
			{
				//printf("YAML_BLOCK_MAPPING_START_TOKEN\n");
				
				break;
			}
			case YAML_SCALAR_TOKEN:
			{
				//Check to see what the type is for the current scalar
				//NOTE: Current scalar is always the type of the previous token
				if (token_old.type == YAML_KEY_TOKEN)
				{
					//Note: Keys will be case sensitive, whereas values false, true, etc. are not
					current_key = (char*) token.data.scalar.value;
					//printf("YAML_SCALAR_TOKEN is a KEY\t%s\n",token.data.scalar.value);
					//std::cout << "YAML_SCALAR_TOKEN is a KEY\t" << current_key << std::endl;
					
					if (doc_start == false)
					{
						std::cout << "\t\tYAML_SCALAR_TOKEN is a DOC_NAME\t" << current_key << std::endl;
					}
					else if (block_key == true && sub_block == false && new_block == true)
					{
						read.key_type = DOCUMENT_HEADER;
						current_block = current_key;
						std::cout << "\t\tYAML_SCALAR_TOKEN is a BLOCK\t" << current_block << std::endl;
					}
					else if (block_key == true && sub_block == true &&  new_block == true)
					{
						read.key_type = SUB_HEADER;
						current_block = current_key;
						std::cout << "\t\tYAML_SCALAR_TOKEN is a SUB-BLOCK\t" << current_block << std::endl;
					}
					else
					{
						read.key_type = KEY;
						std::cout << "\t\tYAML_SCALAR_TOKEN is a KEY\t" << current_key << " and we are in ";
						if (sub_block == false && new_block == false)
						{
							std::cout << "a DOCUMENT\n";
						}
						else if (block_key == false && sub_block == false && new_block == true)
						{
							std::cout << "a HEADER\n";
						}
						else if (block_key == false && sub_block == true && new_block == true)
						{
							std::cout << "a SUB-HEADER\n";
						}
						else
						{
							std::cout << "UNKNOWN\t" << block_key << sub_block << new_block << std::endl;
						}
					}
						
				}
				else if (token_old.type == YAML_VALUE_TOKEN && doc_start == true)
				{
					//First, convert everything to lower case except when case sensitive 
					if (current_block != "masterspecies")
					{
						int i=0;
						while (token.data.scalar.value[i])
						{
							token.data.scalar.value[i] = tolower(token.data.scalar.value[i]);
							i++;
						}
					}

					
					//Must type cast from unsigned char * to char *
					current_value = (char*) token.data.scalar.value;
					d_value = atof(current_value.c_str());
					//printf("YAML_SCALAR_TOKEN is a VALUE\t%s\n",token.data.scalar.value);
					std::cout << "\t\tYAML_SCALAR_TOKEN is a VALUE\t" << current_value << " and we are in ";
					if (sub_block == false && new_block == false)
					{
						std::cout << "a DOCUMENT\n";
					}
					else if (block_key == false && sub_block == false && new_block == true)
					{
						std::cout << "a HEADER\n";
					}
					else if (block_key == false && sub_block == true && new_block == true)
					{
						std::cout << "a SUB-HEADER\n";
					}
					else
					{
						std::cout << "UNKNOWN\t" << block_key << sub_block << new_block << std::endl;
					}
					
					//Check the mapped key for data specifics
					if (current_key == "numvar" || current_key == "num_ssr" || current_key == "num_mbe" || current_key == "num_usr")
						std::cout << "\t\tThis info is for an int!\n";
					if (current_value == "false" || current_value == "true")
						std::cout << "\t\tThis info is for a bool!\n";

					if (d_value == 0)
					{
						if (current_value == "0" || current_value == "0.0")
						{
							std::cout << "\t\tNUMBER\t" << d_value << "\n";
						}
						else
						{
							std::cout << "\t\tNAN\t" << current_value << "\n";
						}
						
					}
					else
					{
						std::cout << "\t\tNUMBER\t" << d_value << "\n";
					}
				}
				else
				{
					std::cout << "\t\t-------ERROR: Unexpected Token------\t" << d_value << "\n";
				}
				
				break;
			}
			case YAML_NO_TOKEN:
			{
				printf("YAML_NO_TOKEN\n");
				return -1; //This is a major error
			}
			case YAML_VERSION_DIRECTIVE_TOKEN:
			{
				//printf("YAML_VERSION_DIRECTIVE_TOKEN\n");
				break;
			}
			case YAML_TAG_DIRECTIVE_TOKEN:
			{
				//printf("YAML_TAG_DIRECTIVE_TOKEN\n");
				break;
			}
			case YAML_DOCUMENT_START_TOKEN:
			{
				printf("YAML_DOCUMENT_START_TOKEN\n--------------------\n");
				doc_start = true;
				break;
			}
			case YAML_DOCUMENT_END_TOKEN:
			{
				printf("YAML_DOCUMENT_END_TOKEN\n---------------------\n\n");
				doc_start = false;
				break;
			}
			case YAML_FLOW_SEQUENCE_START_TOKEN:
			{
				//printf("YAML_FLOW_SEQUENCE_START_TOKEN\n");
				break;
			}
			case YAML_FLOW_SEQUENCE_END_TOKEN:
			{
				//printf("YAML_FLOW_SEQUENCE_END_TOKEN\n");
				break;
			}
			case YAML_FLOW_MAPPING_START_TOKEN:
			{
				//printf("YAML_FLOW_MAPPING_START_TOKEN\n");
				break;
			}
			case YAML_FLOW_MAPPING_END_TOKEN:
			{
				//printf("YAML_FLOW_MAPPING_END_TOKEN\n");
				break;
			}
			case YAML_FLOW_ENTRY_TOKEN:
			{
				//printf("YAML_FLOW_ENTRY_TOKEN\n");
				break;
			}
			case YAML_ALIAS_TOKEN:
			{
				printf("\t---YAML_ALIAS_TOKEN = \t");
				if (doc_start == true && sub_block == false && new_block == false)
				{
					std::cout << "THIS IS ERROR!!!\n";
					break;
				}
				current_alias = (char *) token.data.alias.value;
				std::cout << current_alias << " is a ";
				
				if (doc_start == false)
				{
					std::cout << "DOC-ALIAS = " << current_key << std::endl;
				}
				else if (block_key == true && sub_block == false && new_block == true)
				{
					//current_block = current_key;
					std::cout << "HEAD-ALIAS = " << current_block << std::endl;
				}
				else if (block_key == true && sub_block == true &&  new_block == true)
				{
					//current_block = current_key;
					std::cout << "SUB-ALIAS = " << current_block << std::endl;
				}
				break;
			}
			case YAML_ANCHOR_TOKEN:
			{
				printf("\t---YAML_ANCHOR_TOKEN = \t");
				if (doc_start == true && sub_block == false && new_block == false)
				{
					std::cout << "THIS IS ERROR!!!\n";
					break;
				}
				current_anchor = (char *) token.data.anchor.value;
				std::cout << current_anchor << " is a ";
				
				if (doc_start == false)
				{
					std::cout << "DOC-ANCHOR = " << current_key << std::endl;
				}
				else if (block_key == true && sub_block == false && new_block == true)
				{
					//current_block = current_key;
					std::cout << "HEAD-ANCHOR = " << current_block << std::endl;
				}
				else if (block_key == true && sub_block == true &&  new_block == true)
				{
					//current_block = current_key;
					std::cout << "SUB-ANCHOR = " << current_block << std::endl;
				}
				
				break;
			}
			case YAML_TAG_TOKEN:
			{
				//printf("YAML_TAG_TOKEN\n");
				break;
			}
			default:
			{
				printf("UNKNOWN_TYPE: %d\n", token.type);
				break;
			}
				
		}//END SWITCH CASE
		
		//Set old token object to current token
		token_old = token;
		read.key_type_old = read.key_type;
		read.value_type_old = read.value_type;
		
	} while(token.type != YAML_STREAM_END_TOKEN);
	yaml_token_delete(&token);
	yaml_token_delete(&token_old);
	
	// Cleanup - make sure the file gets closed and the data structure is deleted
	yaml_parser_delete(&parser);
	fclose(fh);
	
	return success;

}
