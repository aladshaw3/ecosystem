/*!
 *  \file mesh.h
 *	\brief Mesh Objects and Associated Sub-Objects
 *	\details This kernel allows for creation of mesh objects from constitutient sub-objects
 *			such as Nodes and Elements (Line, Surface, or Volume). Mesh objects can be used 
 *			to establish physical-chemcial simulations in multi-dimensional space.
 *
 *			\warning This kernel is still under active development. Use with caution!
 *
 *
 *  \author Austin Ladshaw
 *	\date 04/09/2018
 *	\copyright This software was designed and built at the Georgia Institute
 *             of Technology by Austin Ladshaw for Post-Doc research in the area
 *             of adsorption and surface science. Copyright (c) 2018, all
 *             rights reserved.
 */

#ifndef MESH_HPP_
#define MESH_HPP_

#include "macaw.h"

/// Enumeration for the list of valid element types
/** The only types that have been defined are for Boudary and Interior elements.*/
typedef enum {BOUNDARY, INTERIOR} element_type;

/// 3D Vector Object
/** This class structure creates a C++ object for a vector in 3D space. Built using the MACAW matrix
	object, this object contains functions and data associated with working with vectors in 3D space.*/
class Vector3D
{
public:
	Vector3D();									///< Default Constructor
	~Vector3D();								///< Default Destructor
	Vector3D(double x, double y, double z);		///< Construction of RVector with each component
	Vector3D(const Vector3D &v);				///< Copy constructor for the vector
	
	double& operator()(int i);					///< Access to reference of a component of the vector
	double operator()(int i) const;				///< Access to a component of the vector
	double norm();								///< Calculation of the 2-Norm of the vector
	double dot_product(const Vector3D &v);		///< Perform the dot product between two vectors
	double angleRAD(Vector3D &v);				///< Returns the angle between the two vectors in radians
	double angleDEG(Vector3D &v);				///< Returns the angle between the two vectors in degrees
	
	void edit(int i, double value);					///< Editing a single value in the vector
	void set_vector(double x, double y, double z);	///< Editing all values in a vector
	
	Vector3D& operator=(const Vector3D& v);			///< Vector assignment
	Vector3D operator+(const Vector3D& v);			///< Vector Addition
	Vector3D operator-(const Vector3D& v);			///< Vector Subtraction
	double operator*(const Vector3D& v);			///< Vector dot product (short hand = this'*v)
	Vector3D operator*(const double a);				///< Vector-scalar multiplication
	Vector3D operator/(const double a);				///< Vector-scalar division
	Vector3D cross_product(const Vector3D &v);		///< Vector cross product
	
protected:
	Matrix<double> vector;							///< Matrix object to store vector data
	
private:
	
};

/// Node object
/** This class structure creates a C++ object for a node in a mesh. The node will have
	an identifying ID number and a sub_type (either BOUNDARY or INTERIOR) to add in identification
	of different aspects the node has in the overall mesh. All nodes are considered to be points
	in 3D space, whether or not the final mesh is 3D. As such, every node will have a coordinate
	vector (x, y, z) associated with it. */
class Node
{
public:
	Node();												///< Default constructor
	~Node();											///< Default destructor
	
	void DisplayInfo();										///< Display the information associated with this node
	void AssignCoordinates(double x, double y, double z);	///< Assign the coordinates for the node
	void AssignX(double x);									///< Assign the x-coordinate for the node
	void AssignY(double y);									///< Assign the y-coordinate for the node
	void AssignZ(double z);									///< Assign the z-coordinate for the node
	void AssignIDnumber(unsigned int i);					///< Assign the ID number for the node
	void AssignSubType(element_type type);					///< Assign the node_type for the node
	
	bool isSameLocation(Node& node);						///< Returns true if nodes are at same location (witin tolerance)
	bool isSameType(Node& node);							///< Returns true if nodes are of same type
	bool isSameNode(Node& node);							///< Returns true if nodes have same ID number (could indicate error)
	
	double distance(Node& node);							///< Returns the distance between two given nodes
	double angle(Node& n1, Node& n2);						///< Returns angle between n1 and n2 with respect to this node (radians)
	double angle_degrees(Node& n1, Node& n2);				///< Returns angle between n1 and n2 with respect to this node (degrees)
	
	double getX();											///< Returns x value of node
	double getY();											///< Returns y value of node
	double getZ();											///< Returns z value of node
	element_type getType();									///< Returns the type of node
	
	Vector3D operator-(const Node& node);					///< Returns vector formed by difference between nodes
	
protected:
	
private:
	Vector3D coordinates;								///< x, y, z location of the node in space
	unsigned int IDnum;									///< Identification number for the node
	element_type SubType;								///< Sub-type for the node
	double distance_tolerance;							///< Tolerance used to determine if two nodes are in same location
	
};

/// LineElement
/** This class structure creates a C++ object for a line. The line is made up of a reference
	to two distinct nodes. Based on those nodes, we can tell whether or not the line is on
	the boundary of a mesh or part of the interior. We will also know the length of the line
	and the midpoint of the line.*/
class LineElement
{
public:
	LineElement();										///< Default Constructor
	~LineElement();										///< Default Destructor
	
	void DisplayInfo();									///< Print out information to the console
	void AssignNodes(Node &n1, Node &n2);				///< Assign nodes for the line segment
	void AssignIDnumber(unsigned int i);				///< Assign the id number for the line segment
	
	void calculateLength();								///< Calculate and store length value
	void findMidpoint();								///< Find and set the midpoint node
	void determineType();								///< Determine the type of line element
	void evaluateProperties();							///< Calls functions for length and midpoint
	
protected:
	
private:
	Node *node1;										///< Pointer to first node
	Node *node2;										///< Pointer to second node
	double length;										///< Length of the line segment
	Node midpoint;										///< Midpoint node for the line segment
	unsigned int IDnum;									///< Identification number for the line element
	element_type SubType;								///< Type of line segment (BOUDNARY or INTERIOR)

};

/// NodeSet Object
/** This class structure creates a C++ object for a set of nodes. The set of nodes will be a
	listing of all nodes in a domain. Each node must have a unique ID number and unique locations
	within the larger structure. Subsets of the full list are used to construct LineElements,
	SurfaceElements, and VolumeElements of the Mesh object. */
class NodeSet
{
public:
	NodeSet();											///< Default constructor
	~NodeSet();											///< Default destructor
	
	void InitializeNodeSet(unsigned int size);								///< Create a list of nodes of particular size
	void RegisterNode(unsigned int idnum, double x, double y, double z, element_type type);		///< Retister node with specific info
	void RegisterNode(Node &node);											///< Register node from existing node
	
	void AddNode(double x, double y, double z, element_type type);			///< Add node to the list (dynamic)
	
	void AssignDistanceTolerances(double tol);		///< Assign all node tolerances to a specific value
	
protected:
	
private:
	std::vector<Node> node_set;
};

/// Test function for MESH kernel
/** This function runs tests on the mesh objects. */
int MESH_TESTS();

#endif /* MESH_HPP_ */
