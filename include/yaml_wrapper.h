//----------------------------------------
//  Created by Austin Ladshaw on 07/29/15
//  Copyright (c) 2015
//	Austin Ladshaw
//	All rights reserved
//----------------------------------------

/*
 
 DISCLAIMER: Niether Austin Ladshaw, nor the Georgia Institute of Technology, is the author or owner of any YAML Library or source code.
 Only the files labeld "yaml_wrapper" were created by Austin Ladshaw for the sole purpose of running and testing the yaml code before
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

#include "yaml.h"
#include "error.h"
#include <map>
#include <string>
#include <iostream>
#include <utility>
#include <stdexcept>

#ifndef YAML_WRAPPER_HPP_
#define YAML_WRAPPER_HPP_

typedef enum data_type
{ STRING, BOOLEAN, DOUBLE, INT, UNKNOWN } data_type;			//Data type for ValueType

typedef enum header_state { ANCHOR, ALIAS, NONE } header_state;	//State of the Header

class ValueTypePair
{
public:
	ValueTypePair();										//Default constructor
	~ValueTypePair();										//Default destructor
	ValueTypePair(const std::pair<std::string,int> &vt);	//Constructor by pair
	ValueTypePair(std::string value, int type);				//Construction by string and int
	ValueTypePair(const ValueTypePair &vt);					//Copy constructor
	
	ValueTypePair& operator=(const ValueTypePair &vt);			//Equals operator overload
	
	void editValue(std::string value);					//Edits value to pair with UNKOWN type
	void editPair(std::string value, int type);			//Creates a paired Value type
	void findType();									//Determines the data type
	void assertType(int type);							//Forces a specific data type
	void DisplayPair();									//Display the pair
	
	std::string getString();							//Returns the value of the pair as a string
	bool getBool();										//Returns the value of the pair as a bool
	double getDouble();									//Returns the value of the pair as a double
	int getInt();										//Returns the value of the pair as an int
	std::string getValue();								//Returns the value of the pair as it was given
	int getType();										//Returns the type of the pair
	std::pair<std::string,int> &getPair();				//Returns the pair
	
protected:
	
private:
	std::pair<std::string,int> Value_Type;				//Holds the Value and Type info
	int type;											//Type of the value
	
};

class KeyValueMap
{
public:
	KeyValueMap();										//Default constructor
	~KeyValueMap();										//Default destructor
	KeyValueMap(const std::map<std::string,std::string> &map);//Construct from a map of strings
	KeyValueMap(std::string key, std::string value);	//Construct one element in the map
	KeyValueMap(const KeyValueMap &map);				//Copy constructor
	
	KeyValueMap& operator=(const KeyValueMap &map);			//Equals overload
	ValueTypePair& operator[](const std::string key);		//Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	//Return the ValueType at the give key
	
	std::map<std::string, ValueTypePair > & getMap();					//Return a reference to the map object
	std::map<std::string, ValueTypePair>::const_iterator end() const;	//Returns a const iterator pointing to the end of the list
	std::map<std::string, ValueTypePair>::iterator end();				//Returns an iterator pointing to the end of the list
	std::map<std::string, ValueTypePair>::const_iterator begin() const;	//Returns a const iterator pointing to the beginning of the list
	std::map<std::string, ValueTypePair>::iterator begin();				//Returns an iterator pointing to the beginning of the list
	
	void clear();														//Clears the map
	void addKey(std::string key);										//Adds a key to the object with a default value
	void editValue4Key(std::string val, std::string key);				//Edits a given value for a pre-existing key
	void editValue4Key(std::string val, int type, std::string key);		//Edits a value for a pre-existing key and asserts type
	
	void addPair(std::string key, ValueTypePair val);	//Adds a pair object to the map
	void addPair(std::string key, std::string val);		//Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int type); //Adds a pair object and asserts a type
	void findType(std::string key);						//Find what data type the value at the key is
	void assertType(std::string key, int type);			//Assert the given type at the given key
	void findAllTypes();								//Find all types for all data in map
	void DisplayMap();									//Print out the map to console
	
	int size();											//Returns the size of the map
	std::string getString(std::string key);				//Retrieve the string at the key
	bool getBool(std::string key);						//Retrieve the boolean at the key
	double getDouble(std::string key);					//Retrieve the double at the key
	int getInt(std::string key);						//Retrieve the int at the key
	std::string getValue(std::string key);				//Retrieve the value at the key
	int getType(std::string key);						//Retrieve the type at the key
	ValueTypePair& getPair(std::string key);			//Retrieve the pair at the key
	
protected:
	
private:
	std::map<std::string, ValueTypePair > Key_Value;		//Map of Keys and Values paired with types	
};

class SubHeader
{
public:
	SubHeader();										//Default Constructor
	~SubHeader();										//Default Destructor
	SubHeader(const SubHeader &subheader);				//Copy constructor
	SubHeader(const KeyValueMap &map);					//Construction by existing map
	SubHeader(std::string name);						//Construction by name only
	SubHeader(std::string name, const KeyValueMap &map);//Construction by name and map
	
	SubHeader& operator=(const SubHeader &sub);			//Equals overload
	ValueTypePair& operator[](const std::string key);		//Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	//Return the ValueType at the give key
	
	KeyValueMap& getMap();									//Returns reference to the KeyValueMap object
	
	void clear();																//Empty out data contents
	void addPair(std::string key, std::string val);								//Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int type);					//Adds a pair object and asserts a type
	void setName(std::string name);												//Sets the name of the subheader
	void setAlias(std::string alias);											//Set the alias without type specification
	void setAlias(std::string alias, int state);								//Sets the alias of the subheader
	void setNameAliasPair(std::string name, std::string alias, int state);		//Sets the name and alias of the subheader
	void setState(int state);													//Sets the state of the subheader
	void DisplayContents();														//Display the contents of the subheader
	
	std::string getName();								//Return the name of the subheader
	std::string getAlias();								//Return the alias of the subheader, if one exists
	bool isAlias();										//Returns true if subheader is an alias
	bool isAnchor();									//Returns true if subheader is an anchor
	int getState();										//Returns the state of the subheader
	
protected:
	KeyValueMap Data_Map;								//A Map of Keys and Values
	std::string name;									//Name of the subheader
	std::string alias;									//Name of the alias for the subheader
	int state;											//State of the header
	
private:

	
};

class Header : SubHeader
{
public:
	Header();											//Default Constructor
	~Header();											//Default Destructor
	Header(const Header &head);							//Copy constructor
	Header(std::string name);							//Constructor by header name
	Header(const KeyValueMap &map);						//Constructor by existing map
	Header(std::string name, const KeyValueMap &map);	//Constructor by name and map
	Header(std::string key, const SubHeader &sub);		//Constructor by single subheader object
	
	Header& operator=(const Header &head);				//Equals overload
	ValueTypePair& operator[](const std::string key);	//Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	//Return the ValueType at the given key
	SubHeader& operator()(const std::string key);		    //Return the SubHeader reference at the given key
	SubHeader operator()(const std::string key) const;		//Return the SubHeader at the given key
	
	std::map<std::string, SubHeader> & getSubMap();			//Return the reference to the SubHeader Map
	KeyValueMap & getDataMap();								//Return the reference to the KeyValueMap
	SubHeader & getSubHeader(std::string key);				//Return the subheader at the given key
	std::map<std::string, SubHeader>::const_iterator end() const;	//Returns a const iterator pointing to the end of the list
	std::map<std::string, SubHeader>::iterator end();				//Returns an iterator pointing to the end of the list
	std::map<std::string, SubHeader>::const_iterator begin() const;	//Returns a const iterator pointing to the begining of the list
	std::map<std::string, SubHeader>::iterator begin();				//Returns an iterator pointing to the begining of the list
	
	void clear();											//Clear out the SubMap, KeyValueMap, and other info
	void resetKeys();										//Reset the keys of the SubMap
	void changeKey(std::string oldKey, std::string newKey);	//Change one of the keys in the map
	void addPair(std::string key, std::string val);			//Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int t);	//Adds a pair object and asserts a type
	void setName(std::string name);							//Set the name of the Header
	void setAlias(std::string alias);						//Set the alias of the header, if any
	void setNameAliasPair(std::string n, std::string a, int s);	//Set the name, alias, and state for the header
	void setState(int state);								//Set the state of the header, if any
	void DisplayContents();									//Display the contents of the header object
	
	void addSubKey(std::string key);						//Adds a key to the SubHeader Map
	void copyAnchor2Alias(std::string alias, SubHeader &ref);//Find the anchor in the map, and copy to the Header reference given
	
	int size();												//Return the size of the Sub_Map
	std::string getName();									//Return the name of the header
	std::string getAlias();									//Return the alias of the header
	int getState();											//Return the state of the header
	bool isAlias();											//Returns true if the header is an alias
	bool isAnchor();										//Returns true if the header is an anchor
	SubHeader& getAnchoredSub(std::string alias);			//Returns reference to the anchored subheader, if any
	
protected:
	
private:
	std::map<std::string, SubHeader> Sub_Map;			//Map of the contained subheaders in the main header
	
};

class Document : SubHeader
{
public:
	Document();											//Default constructor
	~Document();										//Default destructor
	Document(const Document &doc);						//Copy constructor
	Document(std::string name);							//Constructor by name
	Document(const KeyValueMap &map);					//Constructor by existing map
	Document(std::string name, const KeyValueMap &map);	//Constructor by name and map
	Document(std::string key, const Header &head);		//Constructor by single header
	
	Document& operator=(const Document &doc);				//Equals overload
	ValueTypePair& operator[](const std::string key);		//Return the ValueType reference at the given key
	ValueTypePair operator[](const std::string key) const;	//Return the ValueType at the given key
	Header& operator()(const std::string key);				//Return the Header reference at the given key
	Header operator()(const std::string key) const;			//Return the Header at the given key
	
	std::map<std::string, Header> & getHeadMap();			//Return the reference to the Header Map
	KeyValueMap & getDataMap();								//Return the reference to the KeyValueMap
	Header & getHeader(std::string key);					//Return reference to the Header in map at the key
	std::map<std::string, Header>::const_iterator end() const;	//Returns a const iterator pointing to the end of the list
	std::map<std::string, Header>::iterator end();				//Returns an iterator pointing to the end of the list
	std::map<std::string, Header>::const_iterator begin() const;//Returns a const iterator pointing to the begining of the list
	std::map<std::string, Header>::iterator begin();			//Returns an iterator pointing to the begining of the list
	
	void clear();											//Clear out info in the Document
	void resetKeys();										//Set all keys in the map to match names of the headers
	void changeKey(std::string oldKey, std::string newKey);	//Change a given oldKey in the map to the newKey given
	void revalidateAllKeys();								//Resets and validates keys in header and subheader maps
	void addPair(std::string key, std::string val);			//Adds a pair object to the map (with only strings)
	void addPair(std::string key, std::string val, int t);	//Adds a pair object and asserts a type
	void setName(std::string name);							//Set the name of the Document
	void setAlias(std::string alias);						//Set the alias of the Document
	void setNameAliasPair(std::string n, std::string a, int s);//Set the name, alias, and state of the document
	void setState(int state);								//Set the state of the Document
	void DisplayContents();									//Display the contents of the Document
	
	void addHeadKey(std::string key);						//Add a key to the Header without a header object
	void copyAnchor2Alias(std::string alias, Header &ref);	//Find the anchor in the map, and copy to the Header reference given
	
	int size();												//Return the size of the header map
	std::string getName();									//Return the name of the document
	std::string getAlias();									//Return the alias of the document
	int getState();											//Return the state of the document
	bool isAlias();											//Returns true if the document is an alias
	bool isAnchor();										//Returns true if the document is an anchor
	Header& getAnchoredHeader(std::string alias);			//Returns reference to the anchored header, if any
	Header& getHeadFromSubAlias(std::string alias);			//Returns reference to the Header that contains a Sub with the given alias
	
protected:
	
private:
	std::map<std::string, Header> Head_Map;				//Map of headers contained within the document
	
};

class YamlWrapper
{
public:
	YamlWrapper();										//Default constructor
	~YamlWrapper();										//Default destructor
	YamlWrapper(const YamlWrapper &yaml);				//Copy constructor
	YamlWrapper(std::string key, const Document &doc);	//Constructor by a single document
	
	YamlWrapper& operator=(const YamlWrapper &yaml);	//Equals overload
	Document& operator()(const std::string key);		//Return the Document reference at the given key
	Document operator()(const std::string key) const;	//Return the Document at the given key
	
	std::map<std::string, Document> & getDocMap();		//Return reference to the document map
	Document& getDocument(std::string key);				//Return reference to the document at the key
	std::map<std::string, Document>::const_iterator end() const;	//Returns a const iterator pointing to the end of the list
	std::map<std::string, Document>::iterator end();				//Returns an iterator pointing to the end of the list
	std::map<std::string, Document>::const_iterator begin() const;	//Returns a const iterator pointing to the begining of the list
	std::map<std::string, Document>::iterator begin();				//Returns an iterator pointing to the begining of the list
	
	void clear();										//Clear out the yaml object
	void resetKeys();									//Resets all the keys in DocumentMap to match document names
	void changeKey(std::string oldKey, std::string newKey);	//Change a given oldKey in the map to the newKey given
	void revalidateAllKeys();							//Resets and validates all keys in the structure
	void DisplayContents();								//Display the contents of the wrapper
	
	void addDocKey(std::string key);					//Add a key to the document map
	void copyAnchor2Alias(std::string alias, Document &ref);//Find the anchor in the map, and copy to the Document reference given
	
	int size();											//Return the size of the document map
	Document& getAnchoredDoc(std::string alias);		//Return the reference to the document that is anchored with the given alias
	Document& getDocFromHeadAlias(std::string alias);	//Return reference to the document that contains the header with the given alias
	Document& getDocFromSubAlias(std::string alias);	//Return reference to the document that contains the subheader with the given alias
	
protected:
	
private:
	std::map<std::string, Document> Doc_Map;			//Map of the documents contained within the wrapper
	
};

class yaml_cpp_class
{
public:
	yaml_cpp_class();									//Default constructor
	~yaml_cpp_class();									//Default destructor
	
	int setInputFile(const char *file);					//Set the input file to be read
	int readInputFile();								//Reads through input file and stores into YamlWrapper
	int cleanup();										//Deletes yaml_c objects and closes the input file
	int executeYamlRead(const char *file);				//Runs the full execution of initialization, reading, and cleaning
	YamlWrapper& getYamlWrapper();						//Returns reference to the YamlWrapper Object
	
	void DisplayContents();								//Print out the contents of the read to the console window
	
protected:
	
private:
	YamlWrapper yaml_wrapper;
	FILE *input_file;
	const char *file_name;
	yaml_parser_t token_parser;
	yaml_token_t current_token;
	yaml_token_t previous_token;
	
};

int YAML_WRAPPER_TESTS();

int YAML_CPP_TEST(const char *file);

#endif
